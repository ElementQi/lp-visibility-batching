#include <algorithm>
#include <fstream>
#include <interfaces/highs_c_api.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
int main(int argc, char **argv) {
    if (argc != 3) {
        return 2;
    }
    int m = std::stoi(argv[2]), d = 0;
    std::ifstream input(argv[1]);
    std::vector<int> y;
    std::vector<std::vector<double>> samples;
    std::string line;
    while ((int)y.size() < m && std::getline(input, line)) {
        std::istringstream f(line);
        int label;
        f >> label;
        y.push_back(label);
        std::vector<double> row;
        std::string token;
        while (f >> token) {
            auto pos = token.find(':');
            int j = std::stoi(token.substr(0, pos)) - 1;
            row.resize(std::max<int>(row.size(), j + 1));
            row[j] = std::stod(token.substr(pos + 1));
        }
        d = std::max<int>(d, row.size());
        samples.push_back(row);
    }
    if ((int)y.size() != m) {
        return 2;
    }
    std::vector<HighsInt> starts, indices;
    std::vector<double> v, lo(m, 0), hi(m, 1), cost(m, 1), rl(d + 1, -1), ru(d + 1, 1);
    rl[d] = ru[d] = 0;
    for (int i = 0; i < m; ++i) {
        starts.push_back(v.size());
        for (int j = 0; j < (int)samples[i].size(); ++j) {
            if (samples[i][j]) {
                indices.push_back(j);
                v.push_back(y[i] * samples[i][j]);
            }
        }
        indices.push_back(d);
        v.push_back(y[i]);
    }
    void *h = Highs_create();
    Highs_setBoolOptionValue(h, "output_flag", 0);
    Highs_setIntOptionValue(h, "threads", 1);
    if (Highs_passLp(h, m, d + 1, v.size(), kHighsMatrixFormatColwise, kHighsObjSenseMaximize, 0,
                     cost.data(), lo.data(), hi.data(), rl.data(), ru.data(), starts.data(),
                     indices.data(), v.data()) < 0) {
        return 3;
    }
    if (Highs_run(h) < 0 || Highs_getModelStatus(h) != kHighsModelStatusOptimal) {
        return 4;
    }
    std::vector<double> alpha(m), model(d + 1);
    if (Highs_getSolution(h, alpha.data(), nullptr, nullptr, model.data()) != 0) {
        return 5;
    }
    Highs_destroy(h);
    std::cout << std::setprecision(17);
    for (double x : model) {
        std::cout << x << ' ';
    }
    std::cout << '\n';
    for (double x : alpha) {
        std::cout << x << ' ';
    }
    std::cout << '\n';
}
