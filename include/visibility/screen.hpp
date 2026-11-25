// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace visibility {

enum class Label { unresolved, positive, negative };

namespace detail {
using Integer = std::int64_t;
constexpr int fractional_bits = 20;
constexpr Integer unit = 3 * (1 << 14);  // Normalized witness bounds; 1/3 is exact.

struct Interval { Integer lo, hi; };
using Row = std::vector<Interval>;

inline bool tighten(const Row& row, bool negate, std::vector<Integer>& upper) {
    // The selected coefficient is a lower bound on the true signed coefficient.
    // All supported input and difference endpoints have magnitude below 2^25.
    // Saturating here bounds all sums independently of the number of slots.
    constexpr Integer saturation = (Integer{1} << 25) * unit;
    Integer budget = 0;
    for (std::size_t i = 0; i < row.size(); ++i) {
        Integer a = negate ? -row[i].hi : row[i].lo;
        if (a < 0) {
            budget += (-a) * upper[i];
            if (budget >= saturation) return false;
        }
    }
    bool changed = false;
    for (std::size_t i = 0; i < row.size(); ++i) {
        Integer a = negate ? -row[i].hi : row[i].lo;
        if (a > 0 && budget < a * upper[i]) {
            Integer bound = (budget + a - 1) / a;
            if (bound < upper[i]) {
                upper[i] = bound;
                changed = true;
            }
        }
    }
    return changed;
}

inline void propagate(std::vector<Row> rows, std::vector<Integer>& upper, unsigned passes) {
    const auto n = upper.size();
    // A small projected row is a cheap reference for cancelling shared offsets.
    auto norm = [](const Row& row) {
        Integer peak = 0;
        for (const auto& x : row) peak = std::max({peak, -x.lo, x.hi});
        return peak;
    };
    std::size_t anchor = 0;
    for (std::size_t j = 1; j < rows.size(); ++j)
        if (norm(rows[j]) < norm(rows[anchor])) anchor = j;
    const std::size_t count = rows.size();
    for (std::size_t j = 0; j < count; ++j) {
        if (j == anchor) continue;
        Row difference(n);
        for (std::size_t i = 0; i < n; ++i)
            difference[i] = {rows[j][i].lo - rows[anchor][i].hi,
                             rows[j][i].hi - rows[anchor][i].lo};
        rows.push_back(std::move(difference));
    }
    for (unsigned pass = 0; pass < passes; ++pass) {
        bool changed = false;
        for (const auto& row : rows) {
            changed |= tighten(row, false, upper);
            changed |= tighten(row, true, upper);
        }
        if (!changed) break;
    }
}
// Floating-point elimination only proposes row combinations. A second,
// integer pass re-evaluates the proposed coefficients on the input intervals.
inline std::vector<Row> echelon_rows(const std::vector<std::vector<double>>& basis,
                                     const std::vector<Row>& input) {
    auto work = basis;
    const auto m = work.size(), width = work.front().size();
    std::vector<std::vector<double>> weights(m, std::vector<double>(m));
    for (std::size_t i = 0; i < m; ++i) weights[i][i] = 1;
    std::size_t rank = 0;
    for (std::size_t j = 0; j < width && rank < m; ++j) {
        std::size_t pivot = rank;
        for (std::size_t i = rank+1; i < m; ++i)
            if (std::abs(work[i][j]) > std::abs(work[pivot][j])) pivot = i;
        if (std::abs(work[pivot][j]) < 1e-8) continue;
        std::swap(work[pivot], work[rank]);
        std::swap(weights[pivot], weights[rank]);
        for (std::size_t i = 0; i < m; ++i) {
            if (i == rank) continue;
            double ratio = work[i][j] / work[rank][j];
            if (ratio == 0) continue;
            for (std::size_t k = j+1; k < width; ++k)
                work[i][k] -= ratio * work[rank][k];
            work[i][j] = 0;
            for (std::size_t k = 0; k < m; ++k)
                weights[i][k] -= ratio * weights[rank][k];
        }
        ++rank;
    }
    constexpr Integer denominator = Integer{1} << 20;
    std::vector<Row> result;
    for (std::size_t i = 0; i < rank; ++i) {
        double norm = 0;
        for (double w : weights[i]) norm += std::abs(w);
        if (!std::isfinite(norm) || norm == 0) continue;
        Row row(width, {0, 0});
        Integer remaining = denominator;
        for (std::size_t k = 0; k < m; ++k) {
            double scaled = (weights[i][k]/norm)*denominator;
            Integer q = std::clamp(static_cast<Integer>(scaled), -remaining, remaining);
            if (q == 0) continue;
            remaining -= std::abs(q);
            for (std::size_t j = 0; j < width; ++j) {
                const auto x = input[k][j];
                row[j].lo += q * (q > 0 ? x.lo : x.hi);
                row[j].hi += q * (q > 0 ? x.hi : x.lo);
            }
        }
        for (auto& x : row) {
            x.lo = x.lo/denominator - (x.lo%denominator < 0);
            x.hi = x.hi/denominator + (x.hi%denominator > 0);
        }
        result.push_back(std::move(row));
    }
    return result;
}
}  // namespace detail

// Each supplied binary64 entry must enclose its true paired-basis entry within
// entry_error (including input conversion error). Slots are ordered [+ | -].
// A resolved label excludes visibility on the opposite side. Unresolved labels
// can be passed to the capped LP with the full original feasible box.
inline std::vector<Label> screen(const std::vector<std::vector<double>>& basis,
                                 double entry_error, unsigned passes = 4) {
    if (basis.empty() || basis.front().empty() || basis.front().size() % 2)
        throw std::invalid_argument("expected a nonempty paired basis");
    if (!std::isfinite(entry_error) || entry_error < 0)
        throw std::invalid_argument("expected a finite nonnegative error bound");
    const auto width = basis.front().size(), n = width / 2;
    std::vector<Label> result(n, Label::unresolved);
    bool supported = entry_error <= 1;
    for (const auto& row : basis) {
        if (row.size() != width) throw std::invalid_argument("inconsistent basis width");
        for (double x : row) {
            if (!std::isfinite(x)) throw std::invalid_argument("nonfinite basis entry");
            supported &= std::abs(x) <= 4;
        }
    }
    if (!supported || !passes) return result;
    using namespace detail;
    Integer radius = static_cast<Integer>(std::ceil(std::ldexp(entry_error, fractional_bits)));
    std::vector<Row> input;
    for (const auto& source : basis) {
        Row row;
        row.reserve(width);
        for (double value : source) {
            double x = std::ldexp(value, fractional_bits);
            row.push_back({static_cast<Integer>(std::floor(x)) - radius,
                           static_cast<Integer>(std::ceil(x)) + radius});
        }
        input.push_back(std::move(row));
    }
    std::vector<Integer> upper[2] = {std::vector<Integer>(n, unit), std::vector<Integer>(n, unit)};
    auto screen_rows = [&](const std::vector<Row>& candidates) {
        if (candidates.empty()) return;
        for (int block = 0; block < 2; ++block) {
            std::vector<Row> rows;
            rows.reserve(candidates.size());
            for (const auto& row : candidates) {
                auto first = row.begin()+block*n;
                rows.emplace_back(first, first+n);
            }
            propagate(std::move(rows), upper[block], passes);
            for (std::size_t i = 0; i < n; ++i) {
                if (upper[block][i] < unit/3 && result[i] == Label::unresolved)
                    result[i] = block == 0 ? Label::positive : Label::negative;
            }
        }
    };
    screen_rows(input);
    const auto left = static_cast<std::size_t>(std::count(result.begin(), result.end(),
                                                          Label::unresolved));
    if (left && left <= n/4) screen_rows(echelon_rows(basis, input));
    return result;
}

}  // namespace visibility
