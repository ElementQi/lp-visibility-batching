# Deterministic graph calibration by balanced row elimination

Let $C$ have rank $r$ on $N$ coordinates, with positive price action $C D$ and the reciprocal action on $C^\perp$. A pricing oracle must label each slot visible in a nonnegative vector of either side at threshold $1/3$. The input primitive is an exact minor sign plus a lazy short log proxy, with the common-lift guarantee for the signed cofactor lists at this price. All logarithms in bit budgets are base two. A request can restart the upstream compressor. The reductions count every such request, including requests made before an interruption.

| Construction | Absolute working precision $E$ | Minor lookups per price attempt | Arithmetic on returned short numbers |
|---|---:|---:|---|
| Ordered triangular output, tightened for LP | $O(N+\log N)$ | $O(rN)$ for a one-side flag, $O(N^2)$ for the published paired calibration | Elimination and decoding |
| Full-chart determinant exchanges | $O(\log N)$ | $O(r^2N\log(r+1))$ | Cofactor scans |
| Balanced, refreshed after each merge | $O(\log N)$ | $O(rN(1+\log r))$ | Deterministic basis exchanges on short explicit matrices |
| Balanced, only leaf queries | $O(\log N(1+\log r))$ | $O(rN)$ | Same short-matrix operation |

The last two rows apply to every LP price, including every restriction in witness recovery. The signed-cycle bound gives $K=E+O(\log(N+1)+\log(E+1))$ for the upstream type exponent. Both rows yield an improved deterministic lookup budget even against the linear-precision ordered-array estimate. Recursive pivot choices change the cofactor lists seen at the leaves. This extra freedom explains how the leaf-only procedure can pass the lower bound for decoding a fixed triangular array.

## The exact flag and the merge

Fix any full-rank ordered row basis $A$. In the proof use the unknown effective diagonal for the current price. If $T$ has $a$ pinned columns and the leading $a$ rows on $T$ have nonzero determinant, consider the span of the first $a+k$ priced rows restricted to vectors vanishing on $T$, then omit those coordinates. This space has dimension $k$. Its maximal minors on a $k$-list $J$ are proportional, with a common nonzero factor, to the parent prefix minors on $(T,J)$. Eliminate the leading $a$ rows using their pivots and expand a determinant along $T$ to see the proportionality. At fixed $(a,k,T)$ the signs needed for differences and normalized rows come from these same determinants. A possible factor involving prices on $T$ cancels. Every elementary query in the recursion is therefore a cofactor list of the original row flag.

At $k=1$ ask one list with prefix $T$, choose a proxy winner, and normalize the signed row to pivot one. Its true free entries have magnitude at most $1+8e$, with entry error at most $9e$ after dyadic rounding, where $e=2^{-E}$. At $k>1$ take $l=\lfloor k/2\rfloor$. Recurse on the leading $l$ rows of this space, obtaining pivot columns $A_1$ and a bounded chart. Recurse on the remaining flag with $T$ extended by $A_1$, obtaining pivots $A_2$. In the surviving coordinate order $(A_1,A_2,J)$ the two true row blocks have the shapes

$$
\begin{pmatrix}I&X&Y\\0&I&Z\end{pmatrix}.
$$

Their row spaces are transverse because the first pivot block is identity and the second space vanishes there. Their dimensions sum to $k$, so they span exactly the desired parent. The graph on the concatenated pivots is given by the two free blocks $Y-XZ$ and $Z$. If all three subblocks and their approximations have entries bounded by three and the child errors are at most $h\le1$, the merged entry error is at most $(1+7k)h$. It suffices here that the true entries are bounded by three and the approximations are within $h$. The explicit approximate matrix $G$ has identity pivot columns and all entries $O(k+1)$.

## Stabilize the short matrix before taking another original minor

Apply any single exchange with $|M_{ij}|>2$ on this *explicit short rational* graph until no such entry remains. It doubles the absolute determinant measured in the initial $G$. Hadamard bounds every volume by $(c(k+1)\sqrt{k})^k$ for an absolute $c$. Thus $O(k\log(k+1))$ exchanges suffice. The exchanges use only rational tableau updates. On termination $H=G_B^{-1}G$ has all entries bounded by two. Since the old identity columns occur in $G$, **every entry of $G_B^{-1}$ is also a column entry of $H$**.

Write the true preliminary graph as $G+\Delta$ with $\|\Delta\|_{\max}\le\delta$. Then $F=G_B^{-1}\Delta$ satisfies $\|F\|_{\max}\le2k\delta$. If $2k^2\delta\le1/2$, the row-sum norm of $(I+F_B)^{-1}$ is at most two. In particular the chosen columns are a genuine basis of the true parent. Its true final chart satisfies

$$
H_\star-H=(I+F_B)^{-1}(F-F_BH),\qquad
\|H_\star-H\|_{\max}\le4k(1+2k)\delta.
$$

These estimates use the identity block of the preliminary graph. A generic well-conditioned approximate determinant is unnecessary. The conclusion holds for every common lift compatible with the transcript. This matters when a later source reply further constrains earlier hidden prices. If a lookup restarts, discard this unfinished construction and apply the ordinary restart rule.

## Two error budgets

In the **refreshed** algorithm set $E_0=\lceil\log_2(10^8(N+1)^3)\rceil$. Each child yields $h\le160e$. Substituting $\delta\le(1+7k)h$ into the displayed transfer inequalities gives $2k^2\delta<1/2$ and final error below $1/4$. Hence the true entries on the selected basis are bounded by $9/4$. For every selected pivot request its fundamental cofactor list with all other selected pivots and $T$ fixed. The true pivot normalized against the proxy winner has magnitude at least $1/(9/4)$, since that winner is another entry of this stable true row (or the pivot itself). The computed normalized pivot has magnitude at least $1/4$. The normalized entry errors are $\le9e$. Division and a final dyadic rounding now give entry errors at most $4(1+3)9e+e\le160e$. This closes the induction with a fresh precision budget independent of depth.

In the **leaf-only** algorithm set $B_0=128(N+1)^3$, $L=\lceil\log_2(\max(1,r))\rceil$ and $E=E_0+L\lceil\log_2 B_0\rceil$. Round the stabilized $H$ to the $e$ grid and make zero source requests at the internal node. For $h\ge9e$, the transfer bound and final rounding give an error at most $B_0h$. At height $j$ the error is therefore at most $9e B_0^j\le9\cdot2^{-E_0}$. This bound also ensures all the induction conditions and true entry bounds. The implementation returns the propagated error at the actual root. Both budgets give $6N^2d<1/(64N)$ for the root's error $d$. The capped graph LP therefore returns every required label. A third choice can refresh only at selected depths and charge the product of the amplification factors between consecutive refreshes.

## Requests, bit lengths and a complete-run comparison

Each leaf asks at most $N$ prefix minors. At an internal node of rank $k$ a refresh asks at most $kN$ more. The sum of ranks on any depth is at most $r$, and the maximum depth is $L$. Thus the source counts are at most $rN(1+L)$ with refresh and $rN$ without. Repeated original minors can be cached. Candidate selection may inspect many entries, but each of those operations touches the constructed short array. The exchange kernel uses $O(kN)$ scalar operations per swap, so summing over the balanced tree gives $O(Nr^2\log(r+1))$ short arithmetic operations, in addition to the source, normalization and output. Indeed the sum of squared ranks drops geometrically with depth up to rounding. Product blocks cost $O(k^2N)$ per node. Within a node all iterates are determinant ratios of its fixed preliminary matrix. Dyadic rounding at every return bounds input entry bit lengths at a node by $O(E+\log N)$ (including the product stage); hence intermediate numerator and denominator bit lengths are polynomial by Hadamard. Source log proxies and the price enter only differences, comparisons and bounded exponentiations as in the original calibration.

Hold the published dissociation cap, market, tree depth, nonroot branching and expert count recipe fixed. In the [complete-run ledger](calibration.md#an-explicit-whole-learner-probe-improvement), the scheduled query count is proportional, up to common dimension factors, to $K+\log(D+1)$. Against a strengthened triangular baseline with $K=\Theta(N)$ and the usual $O(N^2)$ minor-request charge, the refreshed route reduces the scheduled probe count by order $N/\log N$ and charges $O(N^2\log N)$ requests per probe. The leaf-only route saves order $N/\log^2 N$ in scheduled calls while keeping $O(N^2)$ requests per probe. In both cases the **resulting full-run minor-request upper budget** improves by order $N/\log^2 N$. Against the original $\Theta(N^3)$ type exponent the corresponding charged improvement is of order $N^3/\log^2 N$. Comparisons of actual work also depend on the costs of each source request and the implicit market computation. Recovering final witnesses can reuse one parent compressor across all restrictions.

The real-matrix exchange idea is classical barycentric-spanner / strong rank-revealing linear algebra. The contribution of this construction is the ordered-prefix oracle reduction, reconditioning on short data, uniform error transfer across compatible lifts, and the deterministic precision/request tradeoff for the strongly polynomial learner. [Gu and Eisenstat](https://doi.org/10.1137/0917055) · [Awerbuch and Kleinberg](https://www.cs.cornell.edu/~rdk/papers/OLSP.pdf) · [Upstream cofactor interface](https://github.com/openai/math/blob/3344da59991ed94f9ba711217b29c1d760badfa6/preprints/Strongly-polynomial-linear-programming-November-24-2026/build/calibration.tex)
