# Faster visibility queries for strongly polynomial linear programming

**[Mengqi Li](https://elementqi.github.io/) · 25 November 2026**

This repository develops complementary graph charts for the visibility queries in OpenAI’s *Strongly polynomial linear programming*. For a general LP price query on $N$ coordinates, a bounded graph chart needs $\Theta(\log(N+1))$ bits of absolute entry precision to distinguish all one-third-visible heads. The ordered triangular array can require $\Theta(N)$ bits, even with every represented entry sign supplied. Both statements have matching upper and lower bounds for their respective input representations. [Proof and separating families](docs/calibration.md)

The chart works with the original online minor compressor. Signed price comparisons form difference constraints on $2N$ vertices. Combining their cycle bound with the new chart reduces the compressor’s affine-type exponent from the published $\Theta(N^3)$ to $O(\log N)$. Holding the published learner and its other budgets fixed, this saves a factor $\Theta(N^3/\log N)$ in the number of scheduled price queries per completed run. Balanced calibration also lowers the full-run minor-request budget against a linear-precision triangular baseline. [Complexity accounting](docs/calibration.md#an-explicit-whole-learner-probe-improvement) · [Upstream algorithm](https://github.com/openai/math/blob/3344da59991ed94f9ba711217b29c1d760badfa6/preprints/Strongly-polynomial-linear-programming-November-24-2026/build/lp.tex)

## Use a graph chart

For an explicit rational row space, build a source once. The chart represents its priced row space in pivot/free order as $[I\; M]$ and its polar as $[-M^\top\; I]$. The same source supplies restricted problems and either side of the alternative.

```python
from visibility.physical import PhysicalSource
from visibility.axis import AxisTable

source = PhysicalSource([[1, 0, -1, 2], [0, 1, 1, 1]])
chart = source.balanced([0, 0, 0, 0])
physical = source.realize(chart)
child = source.chart([0, 0, 0], pinned=[2], polar=True)
axes = AxisTable([0, 0, 0, 0], source.rank, source.view())
far = axes.line(2).chart(10**100)
```

Balanced row elimination uses $O(rN(1+\log r))$ minor requests at logarithmic precision. A leaf-only version uses $O(rN)$ at $O(\log(N+1)(1+\log r))$ bits. When rank is proportional to dimension, a family with opposite mandatory signs attains that precision order for the leaf transcript. [Separating family](docs/leaf-precision.md) An exact rational diagonal realizes the signed intervals on the original input and transports positive-margin solutions back. [Balanced construction](docs/balanced-calibration.md) · [Rational realization](docs/rational-lift.md)

Two independent refinements apply to the general learner. A sharper bound on unqueried proxies saves a factor of the dissociation cap in the scheduled query count. A four-slot parameter choice reduces the probe bound from $O(N^{60}\log(NR+1))$ to $O(N^{48}\log(NR+1))$ at the same demand range $R$. [Extension proof](docs/extension-height.md) · [Learner proof and constants](docs/price-tree.md)

The source evaluates actual rational minors. The general compressed interface accepts an external epoch’s exact signs and rational log proxies. An axis table needs zero new minor requests for arbitrarily distant changes to any one price. [Compressed interface](visibility/calibration.py) · [Axis theorem](docs/price-lines.md)

The C++ graph screen takes a rounded graph block and its entry error. It verifies signs with integer row intervals, then a compact batch LP can handle unresolved coordinates. Against the same graph batch in HiGHS, it won 15 of 20 held-out LP–price pairs (median ratio 2.24×) and 24 of 24 near-optimum Heart SVM–price pairs (3.91×), and 23 of 24 from four other classification datasets (2.35×). Timings start after graph construction. [Implementation](include/visibility/graph.hpp) · [Benchmark, including all regressions](docs/graph-benchmark.md)

Run the exact Python checks with `python3 -m unittest discover -s bench`.

## One batch LP for all visible coordinates

Use the manuscript's paired box for a returned basis. Write $u_j=y_{j+}+\epsilon$ for its nonnegative shifted positive block. Put $\tau=1/(32N)$. Maximize $\sum_j t_j$ subject to the existing box constraints and $0\le t_j\le\min(u_j,\tau)$. Return plus at coordinate $i$ when $u_i\ge\tau/2$.

Suppose a feasible shifted point $v$ has $v_i\ge 1/4$. Given any feasible pair $(u,t)$, set $a=\tau-t_i$. A feasible competitor uses point $(1-a)u+av$ and cap variables $(1-a)t+(a/4)e_i$. Its score gain is $a(1/4-\sum_jt_j)$. Thus, if the objective gap is at most $g$,

$$ (\tau-t_i)(1/4-N\tau)\le g. $$

Every exact maximizer saturates all witnessed heads. At $g\le\tau/32$, each such cap variable still reaches $6\tau/7$. For an opposite witness $w$ with $w_i\ge1/3$, orthogonality to the true positive block and the calibrated entrywise error put $u_i$ strictly below $\tau/2$. Both implications hold for the same LP answer regardless of degeneracy. Signs at invisible heads can differ from the coordinatewise baseline.

## Screen before building the LP

The C++17 filter accepts a calibrated basis with columns ordered $[+\mid-]$ and a uniform absolute error bound per entry. It returns one of three labels per coordinate. If every label is resolved, the caller can skip LP construction. Otherwise keep the full original box and add cap variables only for unresolved coordinates, with the original $\tau$ and gap tolerance. The benchmark implements this path end to end. [Header](include/visibility/screen.hpp) · [Native caller](bench/paired.cpp)

```cpp
#include "visibility/screen.hpp"

auto labels = visibility::screen(paired_basis, entry_error);
```

For a true row $c$ from one block, any normalized nonnegative witness $w$ from the opposite space satisfies $c^\top w=0$. Start with $U_j=1$. If $a_j\le c_j$ and $0\le w_j\le U_j$, then every $a_i>0$ gives

$$
w_i\le \frac{\sum_j \max(-a_j,0)U_j}{a_i}.
$$

The filter rounds this upper bound **up**. Once $U_i<1/3$, the opposite side cannot be visible there, so the current block's sign is valid. Equality stays unresolved by that bound. Both row signs and row differences can tighten several coordinates at once. When at most a quarter of the heads remain, elimination proposes additional row combinations. The filter quantizes their weights and recomputes their intervals in integer arithmetic. Pass counts are fixed, and the optional elimination costs $O(m^2N+m^3)$ arithmetic operations for an $m$-row basis.

Supply an error bound that includes conversion of the basis to binary64. The filter expands entries outward on a dyadic grid. Finite entries outside its supported magnitude of four, or an error bound above one, yield all unresolved. On the residual batch the positive-head proof uses $k\tau\le1/32$ for $k$ unresolved coordinates. The negative-head bound still uses the original paired box. Lean checks the row inequality and the residual cap bound in `interval_row_bound` and `residual_fraction`. The native tests exercise rounding and strict thresholds.

## Check the proof

With [Elan](https://github.com/leanprover/elan) installed, run from this directory:

```sh
lake exe cache get
lake build
```

The project pins Lean 4.33.0 and the mathlib revision in the Lake manifest. The main theorem takes explicit true and approximate matrices, the paper's precision constants, a feasible answer and its objective-gap guarantee. The graph proofs also check signed potentials, cycle-weight transfer, pivot growth, Schur product error and the perturbed-identity solve. Both files report their axioms. [Visibility proofs](BatchedVisibility.lean) · [Graph lemmas](GraphCalibration.lean)

## Reproduce the benchmark

The paired-box benchmark compares the two tests using the same native lp_solve 5.5 library. The coordinate tests reuse a single LP model. Each setting has one warmup and five measured runs with alternating order. The three seed settings repeat the same chain matrix. The 36 settings contain 28 distinct matrices. On full paired 64-dimensional synthetic boxes, the original run gave these pooled medians across 15 observations per row.

| Family | Coordinate tests (ms) | Batch test (ms) | Ratio |
|---|---:|---:|---:|
| Mixed signs | 68.44 | 10.00 | 6.84 |
| Split sign blocks | 50.76 | 9.01 | 5.63 |
| Chain | 50.39 | 8.38 | 6.01 |

The timing includes model assembly, solving and extraction at one price. The benchmarks use exact-product subspaces in mixed bases and a zero lower slack in binary64. The theorem includes the manuscript's nonzero error tolerance. An independent integer simplex calculation checked all 586 required signs among 840 coordinates of the distinct fixtures. Whole-solver performance, including price search and calibration restarts, is open. [Recorded observations](bench/reference.csv)

On Linux with g++, Python 3 and liblpsolve55 installed, adjust the library directory for your installation:

```sh
LP_DIR=/usr/lib/lp_solve
g++ -O2 -std=c++17 bench/paired.cpp -L"$LP_DIR" -Wl,-rpath,"$LP_DIR" -llpsolve55 -o /tmp/lp-visibility-bench
/tmp/lp-visibility-bench > /tmp/lp-visibility-timings.csv
/tmp/lp-visibility-bench --matrices > /tmp/lp-visibility-matrices
python3 bench/verify.py /tmp/lp-visibility-matrices /tmp/lp-visibility-timings.csv
```

## Screened versus full batch

Both paths below receive the same precomputed basis. Timings include screening, model assembly, solving and decoding. Each cell shows the median full-batch time divided by the median screened time at $N=64$ on an Intel Xeon Gold 6240C. The second gauge applies a dense signed row transformation to the same subspaces. [All 360 observations](bench/screen_reference.csv)

| Matrix family | Shared-offset gauge | Varied gauge |
|---|---:|---:|
| Mixed signs | 0.96× | 0.97× |
| Split sign blocks | 22.85× | 4.18× |
| Chain | 16.67× | 3.52× |

At this dimension all split and chain runs resolved all 64 labels before building the LP. The mixed runs resolved zero and paid the additional screening cost. The exact auditor checks both gauges against coordinate optima on the original subspaces. These measurements cover the visibility kernel at one price.

```sh
/tmp/lp-visibility-bench --screen > /tmp/lp-screen.csv
python3 bench/verify_screen.py /tmp/lp-visibility-matrices /tmp/lp-screen.csv
g++ -O2 -std=c++17 bench/screen_test.cpp -o /tmp/lp-screen-test
/tmp/lp-screen-test
```

## Related work and attribution

OpenAI supplies the calibration and global pricing framework. Single-program maximum-support methods, projection/rescaling on complementary subspaces, and optimal supports in truncated cones are relevant earlier work. The result here concerns approximate objective values and two-sided visibility in a fixed, perturbed calibrated box. The research, implementation and Lean proof were developed with AI assistance. [Mehdiloozad et al.](https://arxiv.org/abs/1503.09014) · [Peña and Soheili](https://arxiv.org/abs/2003.08911) · [Ben-Ameur and Maddaloni](https://arxiv.org/abs/2609.37589)

Code and text are available under the [Apache License 2.0](LICENSE).
