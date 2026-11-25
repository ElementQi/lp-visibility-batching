# Fixed-price tests on LP and classification data

Both paths start with the same binary64 matrix $M$ and error bound $2^{-50}$. The positive and polar spaces are $(a,aM)$ and $(-bM^\top,b)$. The compact baseline maximizes the capped objective in HiGHS using only the positive pivot variables and one cap variable per head. The proposed path first checks graph basis rows, then uses projected rows to propose integer interval certificates. It passes any unresolved heads to the same compact LP. A third mode runs basis-row screening alone. Every path times assembly, screening where enabled, solution and decoding. Exact graph construction precedes the timer.

| Suite | Distinct case–price pairs | Screen faster | Median paired speed ratio | Range |
|---|---:|---:|---:|---:|
| Development LP and Heart kernels | 12 | 10 | 5.38× | 0.53–18.57× |
| Held-out public LP kernels | 20 | 15 | 2.24× | 0.24–12.73× |
| Heart, certified budgets around the optimum | 24 | 24 | 3.91× | 1.36–18.53× |
| Cancer, digits, iris and wine, certified budgets | 24 | 23 | 2.35× | 0.89–12.39× |

Each ratio is the median of five baseline timings divided by the median of five screened timings for that case. The headline is the median of those ratios. Five successive repetitions rotate through the three execution orders, including the first repetition. The runs use one HiGHS thread, version 1.12.0, on an Intel Xeon Gold 6240C. Every projected-screen run in this table resolved all heads before building a residual LP. For RGN at zero price the baseline took 9.66 ms and the screened path took 40.19 ms. Distinct prices on a model and nested Heart subsets share input structure. [All case medians, including the row-only mode](../bench/graph_reference.csv)

## Inputs and witnesses

The development models use the COIN-OR sample collection and the LIBSVM Heart data. For the held-out suite we fixed twelve HiGHS model names and price seeds 0 and 19. The slot cap was 550. Refinery (705 slots) and standata (1,363) exceeded it. The ten models timed were adlittle, avgas, bell5, egout, flugpl, israel, lseu, qap04, rgn and woodinfe. We homogenized the linear constraints, represented free columns as signed differences, added slacks, and relaxed integrality. The seeded prices are independent draws in $[-4,4]$ using Python's Random. The input builder uses exact rational arithmetic and pivots until every graph entry has magnitude at most one. Binary64 conversion is covered by the supplied error bound. [COIN-OR snapshot](https://github.com/coin-or-tools/Data-Sample/tree/a75e56df8b3f44b09fcf7102546a3ca3cbf53a4f) · [HiGHS instances and source](https://github.com/ERGO-Code/HiGHS/tree/755a8e027a99a8d4ecf153a8dde4b2a767cdf384) · [LIBSVM Heart](https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/heart_scale)

For classification we use the epigraph decision problem for $\min_{w,b}\|w\|_1+\sum_i[1-y_i(w^\top x_i+b)]_+$. The intercept is free. For Heart we used the first 32, 64, 128 and 270 observations. The other four datasets use 32 records from each selected class, with a fixed shuffle. Within the selected records each nonconstant feature was rescaled to $[-1,1]$ and rounded to a $1/256$ grid. This setup studies the LP used for training. Prediction accuracy plays no role in the visibility contract. [Subset construction](../bench/ml_subset.py) · [scikit-learn bundled datasets](https://scikit-learn.org/stable/datasets/toy_dataset.html)

Each subset has a rational primal upper bound $U$ from a classifier and hinge slacks. A dual vector supplies a lower bound $L$ after rational rounding, balancing the two label masses and enforcing $\|\sum_i\alpha_i y_i x_i\|_\infty\le1$. We require $U-L<10^{-4}$. Budgets are $\lfloor990L\rfloor/1000$ and $\lceil1010U\rceil/1000$. There are three prices on the target coordinate, $-q,0,q$, with $q=2\operatorname{bitlength}(m+1)+4$ for $m$ samples. On the feasible side the classifier extends to a nonnegative homogenized primal vector of target one. On the infeasible side the budget row minus the $\alpha$-weighted margin rows gives a nonnegative polar vector with target $L-\mathrm{budget}>0$. The builder checks each witness against the rational rows and measures its target after reciprocal price scaling. The target is mandatory when its normalized magnitude is at least one third. Across the two budget suites, all 315 mandatory target checks passed (21 distinct cases × three modes × five repetitions). A single witness below threshold leaves the obligation at that slot undecided.

Separately, fraction-free simplex optimized *every* positive and negative head on five exact graphs. AFIRO at two prices, GALENET, P0033 and Heart32 had 304 coordinates, including 181 mandatory signs. Every tested mask respected those signs. The general guarantee for the screen is its integer row-bound proof. HiGHS reports a numerical optimum when the batch is used. [Exact coordinate checker](../bench/verify_graph.py)

## Reproduction

Build the pinned HiGHS source with shared libraries. For example, set HIGHS to the source checkout and HB to its CMake build directory. Run the following from the repository root.

```sh
g++ -O3 -std=c++17 -I"$HIGHS/highs" -I"$HB" bench/graph.cpp \
    -L"$HB/lib" -Wl,-rpath,"$HB/lib" -lhighs -o /tmp/graph-price
g++ -O3 -std=c++17 -I"$HIGHS/highs" -I"$HB" bench/svm_reference.cpp \
    -L"$HB/lib" -Wl,-rpath,"$HB/lib" -lhighs -o /tmp/svm-reference
mkdir -p /tmp/graph-fixtures
python3 bench/graph_inputs.py "$HIGHS/check/instances/adlittle.mps" adlittle_19 \
    /tmp/graph-fixtures/adlittle_19 --price-seed 19
/tmp/graph-price /tmp/graph-fixtures/adlittle_19
```

With the Heart file at HEART, reproduce its two 64-observation budgets and all three prices.

```sh
mkdir -p /tmp/heart-sol /tmp/heart-graphs
/tmp/svm-reference "$HEART" 64 > /tmp/heart-sol/visibility_svm_solution_64
python3 bench/svm_budgets.py "$HEART" /tmp/heart-sol /tmp/heart-graphs 64
/tmp/graph-price /tmp/heart-graphs/heart64_*
```

For the other datasets first generate the LIBSVM-format input using scikit-learn (reference version 1.2.2), then use the same budget builder with the appropriate tag.

```sh
python3 bench/ml_subset.py digits /tmp/digits64.data
/tmp/svm-reference /tmp/digits64.data 64 > /tmp/heart-sol/visibility_svm_solution_64
python3 bench/svm_budgets.py --tag-prefix digits /tmp/digits64.data /tmp/heart-sol /tmp/digits-graphs 64
/tmp/graph-price /tmp/digits-graphs/digits64_*
```
