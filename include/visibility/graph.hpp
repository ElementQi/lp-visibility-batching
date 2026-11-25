// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "screen.hpp"
#include <numeric>

namespace visibility {
namespace graph_detail {
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;
using detail::Integer;
using detail::Interval;
using detail::Row;

inline double dot(const Vector &a, const Vector &b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.);
}

// Orthogonally project onto a graph with basis [I M]. The smaller of the
// primal and polar graphs is used, so only that many normal equations enter.
class Projector {
    Matrix graph_, triangular_;
    bool polar_;

  public:
    explicit Projector(const Matrix &matrix) : polar_(matrix.front().size() < matrix.size()) {
        if (polar_) {
            graph_.assign(matrix.front().size(), Vector(matrix.size()));
            for (size_t i = 0; i < matrix.size(); ++i) {
                for (size_t j = 0; j < matrix[i].size(); ++j) {
                    graph_[j][i] = -matrix[i][j];
                }
            }
        } else {
            graph_ = matrix;
        }
        const size_t r = graph_.size();
        triangular_.assign(r, Vector(r));
        for (size_t i = 0; i < r; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                double v = (i == j ? 1. : 0.) + dot(graph_[i], graph_[j]);
                for (size_t k = 0; k < j; ++k) {
                    v -= triangular_[i][k] * triangular_[j][k];
                }
                triangular_[i][j] = i == j ? std::sqrt(v) : v / triangular_[j][j];
            }
        }
    }
    void split(const Vector &z, Vector &positive, Vector &negative) const {
        size_t r = graph_.size(), d = graph_.front().size();
        size_t piv = polar_ ? d : 0, free = polar_ ? 0 : r;
        auto &p = polar_ ? negative : positive;
        auto &q = polar_ ? positive : negative;
        for (size_t i = 0; i < r; ++i) {
            double v = z[piv + i] +
                       std::inner_product(graph_[i].begin(), graph_[i].end(), z.begin() + free, 0.);
            for (size_t j = 0; j < i; ++j) {
                v -= triangular_[i][j] * p[piv + j];
            }
            p[piv + i] = v / triangular_[i][i];
        }
        for (size_t i = r; i-- > 0;) {
            double v = p[piv + i];
            for (size_t j = i + 1; j < r; ++j) {
                v -= triangular_[j][i] * p[piv + j];
            }
            p[piv + i] = v / triangular_[i][i];
        }
        std::fill(p.begin() + free, p.begin() + free + d, 0.);
        for (size_t i = 0; i < r; ++i) {
            for (size_t j = 0; j < d; ++j) {
                p[free + j] += p[piv + i] * graph_[i][j];
            }
        }
        for (size_t j = 0; j < z.size(); ++j) {
            q[j] = z[j] - p[j];
        }
    }
};

class RowBounds {
    static constexpr Integer grid = Integer{1} << 30;
    static constexpr Integer weights = Integer{1} << 26;
    static constexpr Integer row_grid = Integer{1} << 20;
    const size_t r_, d_, n_;
    std::vector<Row> entries_;
    std::vector<Integer> upper_[2];
    std::vector<Row> bank_[2];
    size_t left_;

    void update_labels(int side) {
        for (size_t i = 0; i < n_; ++i) {
            if (labels[i] == Label::unresolved && upper_[side][i] < detail::unit / 3) {
                labels[i] = side ? Label::negative : Label::positive;
                --left_;
            }
        }
    }
    static Integer floor_div(Integer value, Integer den) { return value / den - (value % den < 0); }

  public:
    std::vector<Label> labels;
    RowBounds(const Matrix &matrix, double error)
        : r_(matrix.size()), d_(matrix.front().size()), n_(r_ + d_), entries_(r_, Row(d_)),
          left_(n_), labels(n_, Label::unresolved) {
        const Integer radius = static_cast<Integer>(std::ceil(std::ldexp(error, 30)));
        for (size_t i = 0; i < r_; ++i) {
            for (size_t j = 0; j < d_; ++j) {
                double x = std::ldexp(matrix[i][j], 30);
                entries_[i][j] = {static_cast<Integer>(std::floor(x)) - radius,
                                  static_cast<Integer>(std::ceil(x)) + radius};
            }
        }
        for (auto &v : upper_) {
            v.assign(n_, detail::unit);
        }
    }
    bool done() const { return left_ == 0; }
    void replay() {
        for (int side = 0; side < 2; ++side) {
            for (const auto &row : bank_[side]) {
                detail::tighten(row, false, upper_[side]);
            }
            update_labels(side);
        }
    }
    void propose(const Vector &row, int side) {
        const size_t rank = side ? d_ : r_, width = side ? r_ : d_;
        const size_t piv = side ? r_ : 0, free = side ? 0 : r_;
        const auto &upper = upper_[side];
        double loss = 0, head = 0, norm = 0;
        for (size_t j = 0; j < n_; ++j) {
            loss += std::max(0., -row[j]) * upper[j];
            if (labels[j] == Label::unresolved) {
                head = std::max(head, row[j] * upper[j]);
            }
        }
        if (!(head > loss)) {
            return;
        }
        for (size_t i = 0; i < rank; ++i) {
            norm += std::abs(row[piv + i]);
        }
        if (!(norm > 0) || !std::isfinite(norm)) {
            return;
        }
        std::vector<Integer> lo(n_);
        Integer remaining = weights;
        for (size_t i = 0; i < rank; ++i) {
            Integer w = std::clamp(static_cast<Integer>(row[piv + i] / norm * weights), -remaining,
                                   remaining);
            remaining -= std::abs(w);
            if (!w) {
                continue;
            }
            lo[piv + i] = w * grid;
            for (size_t j = 0; j < width; ++j) {
                auto interval = side ? entries_[j][i] : entries_[i][j];
                Integer c = side ? -w : w;
                lo[free + j] += c * (c > 0 ? interval.lo : interval.hi);
            }
        }
        Row certified(n_);
        for (size_t i = 0; i < n_; ++i) {
            Integer v = floor_div(lo[i], weights * (grid / row_grid));
            certified[i] = {v, v};
        }
        if (detail::tighten(certified, false, upper_[side])) {
            update_labels(side);
            if (bank_[side].size() < 2 * n_) {
                bank_[side].push_back(std::move(certified));
            }
        }
    }
    void basis_rows() {
        for (int side = 0; side < 2; ++side) {
            size_t rank = side ? d_ : r_, width = side ? r_ : d_;
            size_t piv = side ? r_ : 0, free = side ? 0 : r_;
            std::vector<Row> rows;
            for (size_t i = 0; i < rank; ++i) {
                Row row(n_, {0, 0});
                row[piv + i] = {row_grid, row_grid};
                for (size_t j = 0; j < width; ++j) {
                    auto a = side ? entries_[j][i] : entries_[i][j];
                    Integer l = side ? -a.hi : a.lo, h = side ? -a.lo : a.hi;
                    row[free + j] = {floor_div(l, grid / row_grid),
                                     -floor_div(-h, grid / row_grid)};
                }
                rows.push_back(std::move(row));
            }
            detail::propagate(std::move(rows), upper_[side], 4);
            update_labels(side);
        }
    }
};
} // namespace graph_detail

// Entries describe [I M] and its polar [-M^T I] in the same pivot/free order.
// entry_error bounds every entry of M relative to the true graph. A finite
// search budget is optional; callers can pass unresolved heads to the graph LP.
inline std::vector<Label> screen_graph(const std::vector<std::vector<double>> &matrix,
                                       double entry_error, size_t steps) {
    if (matrix.empty() || matrix.front().empty()) {
        throw std::invalid_argument("expected a nonempty graph block");
    }
    if (!std::isfinite(entry_error) || entry_error < 0) {
        throw std::invalid_argument("invalid graph entry error");
    }
    const size_t d = matrix.front().size(), n = matrix.size() + d;
    bool supported = entry_error <= 1;
    for (const auto &row : matrix) {
        if (row.size() != d) {
            throw std::invalid_argument("inconsistent graph width");
        }
        for (double x : row) {
            if (!std::isfinite(x)) {
                throw std::invalid_argument("nonfinite graph entry");
            }
            supported = supported && std::abs(x) <= 4;
        }
    }
    if (!supported) {
        return std::vector<Label>(n, Label::unresolved);
    }
    graph_detail::RowBounds bounds(matrix, entry_error);
    bounds.basis_rows();
    if (bounds.done() || !steps) {
        return bounds.labels;
    }
    graph_detail::Projector projector(matrix);
    graph_detail::Vector p(n), q(n);
    for (int side = 0; side < 2; ++side) {
        graph_detail::Vector z(n, 1.);
        for (size_t k = 0; k < steps; ++k) {
            projector.split(z, p, q);
            if (k < 4 || k % 4 == 3) {
                bounds.propose(p, 0);
                bounds.propose(q, 1);
                if (k % 16 == 15) {
                    bounds.replay();
                }
                if (bounds.done()) {
                    return bounds.labels;
                }
            }
            const auto &projected = side ? q : p;
            for (size_t i = 0; i < n; ++i) {
                double floor = bounds.labels[i] == Label::unresolved ? 1. : 0.;
                z[i] += std::max(floor, 2 * projected[i] - z[i]) - projected[i];
            }
        }
    }
    return bounds.labels;
}
} // namespace visibility
