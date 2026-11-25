// SPDX-License-Identifier: Apache-2.0
#include "../include/visibility/graph.hpp"
#include <chrono>
#include <fstream>
#include <interfaces/highs_c_api.h>
#include <iostream>
#include <memory>
#include <string>

using V = std::vector<double>;
using Matrix = std::vector<V>;
using Labels = std::vector<int>;
using Clock = std::chrono::steady_clock;
using visibility::Label;

static void require(bool ok) {
    if (!ok) {
        throw std::runtime_error("HiGHS failed to solve the graph batch");
    }
}

Labels graph_lp(const Matrix &m, double error, const std::vector<Label> &signed_labels) {
    const int r = m.size(), d = m.front().size(), n = r + d;
    std::vector<int> active;
    Labels answer(n);
    for (int i = 0; i < n; ++i) {
        if (signed_labels[i] == Label::unresolved) {
            active.push_back(i);
        } else {
            answer[i] = signed_labels[i] == Label::positive;
        }
    }
    if (active.empty()) {
        return answer;
    }
    const int cols = r + active.size(), rows = d + active.size();
    const double pad = n * error, tau = 1. / (32 * n);
    V cost(cols), lo(cols, 0), hi(cols, 1), rl(rows, -pad), ru(rows, 1 + pad), val;
    std::vector<HighsInt> starts(rows), idx;
    for (int j = 0; j < rows; ++j) {
        starts[j] = idx.size();
        int slot = j < d ? r + j : active[j - d];
        if (slot < r) {
            idx.push_back(slot);
            val.push_back(1);
        } else {
            for (int i = 0; i < r; ++i) {
                if (m[i][slot - r]) {
                    idx.push_back(i);
                    val.push_back(m[i][slot - r]);
                }
            }
        }
        if (j >= d) {
            int k = r + j - d;
            hi[k] = tau;
            cost[k] = 1;
            idx.push_back(k);
            val.push_back(-1);
            ru[j] = 1e30;
        }
    }
    std::unique_ptr<void, decltype(&Highs_destroy)> h(Highs_create(), Highs_destroy);
    Highs_setBoolOptionValue(h.get(), "output_flag", 0);
    Highs_setIntOptionValue(h.get(), "threads", 1);
    Highs_setDoubleOptionValue(h.get(), "primal_feasibility_tolerance", 1e-9);
    Highs_setDoubleOptionValue(h.get(), "dual_feasibility_tolerance", 1e-9);
    require(Highs_passLp(h.get(), cols, rows, val.size(), kHighsMatrixFormatRowwise,
                         kHighsObjSenseMaximize, 0, cost.data(), lo.data(), hi.data(), rl.data(),
                         ru.data(), starts.data(), idx.data(), val.data()) >= 0);
    require(Highs_run(h.get()) >= 0 && Highs_getModelStatus(h.get()) == kHighsModelStatusOptimal);
    V x(cols);
    require(Highs_getSolution(h.get(), x.data(), nullptr, nullptr, nullptr) == 0);
    for (size_t j = 0; j < active.size(); ++j) {
        answer[active[j]] = x[r + j] > tau / 2;
    }
    return answer;
}

int main(int argc, char **argv) {
    std::cout << "case,n,rep,mode,screen_ms,total_ms,left,status,labels\n";
    if (argc < 2) {
        std::cerr << "usage: graph CASE [CASE ...]\n";
        return 2;
    }
    bool failed = false;
    for (int z = 1; z < argc; ++z) {
        std::ifstream f(argv[z]);
        std::string name, v;
        int r, d;
        if (!(f >> name >> r >> d) || r < 1 || d < 1) {
            return 2;
        }
        Matrix m(r, V(d));
        for (auto &row : m) {
            for (auto &x : row) {
                if (!(f >> v)) {
                    return 2;
                }
                x = std::stod(v);
            }
        }
        int n = r + d;
        for (int rep = 0; rep < 5; ++rep) {
            for (int j = 0; j < 3; ++j) {
                int mode = (j + rep) % 3;
                const auto start = Clock::now();
                auto labels = mode ? visibility::screen_graph(m, 0x1p-50, mode == 2 ? 8 * n : 0)
                                   : std::vector<Label>(n, Label::unresolved);
                const auto mid = Clock::now();
                Labels out;
                std::string status = "ok";
                try {
                    out = graph_lp(m, 0x1p-50, labels);
                } catch (const std::exception &) {
                    status = "full_retry";
                    try {
                        out = graph_lp(m, 0x1p-50, std::vector<Label>(n, Label::unresolved));
                    } catch (const std::exception &) {
                        status = "failure";
                        failed = true;
                    }
                }
                const auto end = Clock::now();
                auto ms = [](auto a, auto b) {
                    return std::chrono::duration<double, std::milli>(b - a).count();
                };
                std::cout << name << ',' << n << ',' << rep << ',' << mode << ',' << ms(start, mid)
                          << ',' << ms(start, end) << ','
                          << std::count(labels.begin(), labels.end(), Label::unresolved) << ','
                          << status << ',';
                for (auto x : out) {
                    std::cout << x;
                }
                std::cout << std::endl;
            }
        }
    }
    return failed ? 1 : 0;
}
