# Batched visibility for linear programming

**[Mengqi Li](https://elementqi.github.io/) · v0.1.0 · 25 November 2026**

A single linear program can label every required coordinate of a calibrated price query in OpenAI's *Strongly polynomial linear programming*. The published test solves one feasibility problem per coordinate. The batch program adds one variable per coordinate and permits additive objective error. [Original pricing algorithm](https://github.com/openai/math/blob/3344da59991ed94f9ba711217b29c1d760badfa6/preprints/Strongly-polynomial-linear-programming-November-24-2026/build/lp.tex)

## The test

Use the manuscript's paired box for a returned basis. Write $u_j=y_{j+}+\epsilon$ for its nonnegative shifted positive block. Put $\tau=1/(32N)$. Maximize $\sum_j t_j$ subject to the existing box constraints and $0\le t_j\le\min(u_j,\tau)$. Return plus at coordinate $i$ when $u_i\ge\tau/2$.

Suppose a feasible shifted point $v$ has $v_i\ge 1/4$. Given any feasible pair $(u,t)$, set $a=\tau-t_i$. A feasible competitor uses point $(1-a)u+av$ and cap variables $(1-a)t+(a/4)e_i$. Its score gain is $a(1/4-\sum_jt_j)$. Thus, if the objective gap is at most $g$,

$$ (\tau-t_i)(1/4-N\tau)\le g. $$

Every exact maximizer saturates all witnessed heads. At $g\le\tau/32$, each such cap variable still reaches $6\tau/7$. For an opposite witness $w$ with $w_i\ge1/3$, orthogonality to the true positive block and the calibrated entrywise error put $u_i$ strictly below $\tau/2$. Both implications hold for the same LP answer regardless of degeneracy. Signs at invisible heads can differ from the coordinatewise baseline.

## Check the proof

With [Elan](https://github.com/leanprover/elan) installed, run from this directory:

```sh
lake exe cache get
lake build
```

The project pins Lean 4.33.0 and the mathlib revision in the Lake manifest. The main theorem takes explicit true and approximate matrices, the paper's precision constants, a feasible answer and its objective-gap guarantee. The file reports the axioms used by the main results. [Theorem and proof](BatchedVisibility.lean)

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

## Related work and attribution

OpenAI supplies the calibration and global pricing framework. Single-program maximum-support methods, projection/rescaling on complementary subspaces, and optimal supports in truncated cones are relevant earlier work. The result here concerns approximate objective values and two-sided visibility in a fixed, perturbed calibrated box. The research, implementation and Lean proof were developed with AI assistance. [Mehdiloozad et al.](https://arxiv.org/abs/1503.09014) · [Peña and Soheili](https://arxiv.org/abs/2003.08911) · [Ben-Ameur and Maddaloni](https://arxiv.org/abs/2609.37589)

Code and text are available under the [Apache License 2.0](LICENSE).
