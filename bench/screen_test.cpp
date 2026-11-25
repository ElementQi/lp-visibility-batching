// SPDX-License-Identifier: Apache-2.0
#include "../include/visibility/screen.hpp"
#include <array>
#include <cassert>
#include <iostream>
#include <random>

using visibility::Label;
using Matrix = std::vector<std::vector<double>>;

int main() {
    // Head 0 has an opposite witness exactly at 1/3. Ignoring the known
    // perturbation could certify the wrong sign at that threshold.
    Matrix edge{{3, -1.+1./32, 0, 0}, {0, 0, 1, 3}};
    assert(visibility::screen(edge, 0)[0] == Label::positive);
    assert(visibility::screen(edge, 1./32)[0] != Label::positive);
    edge[0][1] = -1;
    assert(visibility::screen(edge, 0)[0] != Label::positive);
    assert(visibility::screen(edge, 2)[0] == Label::unresolved);
    assert(visibility::screen(edge, 0, 0)[0] == Label::unresolved);
    // Arbitrary mixed row gauges and dyadic input error around known lines.
    std::mt19937 rng(301);
    int checked = 0;
    for (int a=-3; a<=3; ++a) for (int b=-3; b<=3; ++b) {
        if (!a && !b) continue;
        for (int rep=0; rep<20; ++rep) {
            Matrix basis(2,std::vector<double>(4));
            for(int k=0;k<2;++k) {
                double x=k==0?0.5:(int(rng()%9)-4)/16.;
                double y=k==1?0.5:(int(rng()%9)-4)/16.;
                basis[k]={x*a,x*b,-y*b,y*a};
                for(auto& v:basis[k]) v+=(int(rng()%3)-1)/1024.;
            }
            const auto labels=visibility::screen(basis,1./1024.);
            for(int block=0;block<2;++block){
                const std::array<int,2> v=block?std::array<int,2>{-b,a}:std::array<int,2>{a,b};
                if(v[0]*v[1]<0) continue;
                const int peak=std::max(std::abs(v[0]),std::abs(v[1]));
                for(int i=0;i<2;++i) if(3*std::abs(v[i])>=peak) {
                    assert(labels[i]!=(block?Label::positive:Label::negative));
                    ++checked;
                }
            }
        }
    }
    // Direct sums of lines in dense, diagonally dominant row gauges.
    // This also exercises the optional row-combination pass at larger n.
    for (int rep = 0; rep < 80; ++rep) {
        constexpr int n = 8;
        Matrix truth(n, std::vector<double>(2*n));
        std::array<Label,n> required;
        required.fill(Label::unresolved);
        for (int i = 0; i < n; i += 2) {
            int b = int(rng()%7)-3;
            truth[i][i]=1; truth[i][i+1]=b;
            truth[i+1][n+i]=-b; truth[i+1][n+i+1]=1;
            if (b >= 0) required[i]=Label::positive;
            if (b > 0) required[i+1]=Label::positive;
            if (b < 0) required[i]=Label::negative;
            if (b <= 0) required[i+1]=Label::negative;
        }
        Matrix observed(n, std::vector<double>(2*n));
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n; ++k) {
                double w=(i==k?.5:0.)+(int(rng()%3)-1)/64.;
                for(int j=0;j<2*n;++j) observed[i][j]+=w*truth[k][j];
            }
            for(double& x:observed[i]) x+=(int(rng()%3)-1)/65536.;
        }
        const auto labels=visibility::screen(observed,1./65536.);
        for(int i=0;i<n;++i) {
            assert(labels[i]==Label::unresolved || labels[i]==required[i]);
            ++checked;
        }
    }
    std::cout<<checked<<" nonzero-error visible heads checked\n";
}
