# The global learner, expert cells and two smaller deterministic probe trees

The bounded graph result changes calibration. The following refinement changes the recursive price learner itself. The four-slot version retains the original market, investment rate and update exactly. It reduces the published probe bound $O(N^{60}\log(NR+1))$ to $O(N^{48}\log(NR+1))$. A five-slot version gives $O(N^{39}\log(NR+1))$ with a larger implicit market. The ratios of these displayed exponents describe upper bounds. Compared at fixed $R$, the actual prescribed completed tree has the power $\log_g r$ in its branching term, where $g$ is the growth factor below. These schedules concern every signed extremal-price instance, independent of how the labels are acquired.

## General constants in the published proof

Use $k$ iid draws, boost the first $\delta$ of sorted mass by $\beta>1$, and select a heavy pivot from the quantile band $[\gamma,\delta)$. Define $b=\beta/(1+(\beta-1)\delta)$ and $g=\sqrt\beta/(1+(\beta-1)\delta)>1$. The least positive odd depth $h$ with $g^h>2N$ ensures that a leaf median meets its threshold. Any violating coordinate would have received the full boost at more than half the levels, giving final mass greater than one. Keep the published share $c$, market rate $\eta$, and bonus $L=\ln(1+(c/2)\mu)$, or substitute the rational values below. The current proof experts have mass one. Reserve an additional mass $u\lambda^{j-l}$ for each future pivot level $j\ge l$, so total mass is at most $M=1+u/(1-\lambda)$. These reserves are proof weights. They do not enter actual wealth updates.

The child's current class costs at most $b^k$ times the parent class on a light probe, and at most $b^k/[u(\delta-\gamma)]$ times the matching reserve in expectation on a heavy probe. A future reserve costs at most $b^{k+1}/\lambda$. Class supports are disjoint by pivot-level tag. Thus one may choose a heavy pivot satisfying the bonus-transfer bound with

$$C=b^k\max\{1,[u(\delta-\gamma)]^{-1},b/\lambda\}.$$

For a successful heavy probe put $p=1-(1-\gamma)^k$. Its expected log-market gain is at least $G=p\ln(1+\eta)+(1-p)\ln(1-\eta)$. On a successful light probe it is nonnegative. As in the original proof a failed node has at least half its probes either heavy or child-failed. Let $\epsilon,\sigma>0$ and odd nonroot arity $r$ satisfy

$$G\ge\epsilon,\qquad \ln(1-\eta)+\sigma r/C\ge\epsilon,$$
$$r(\epsilon/2+M\ln(1-c)-\sigma)\ge M\ln(2/c).$$

The first condition pays for successful heavy probes, the second pays for failed children by induction and transfer, and the third pays the parent wealth split and initial log fee. Each failed node with $T\ge r$ probes then has $\sum\pi L\ge\sigma T$. At the root use $r_0=r+s\lceil\log_2(|E|+1)\rceil$ for even $s$ with $s\sigma\ge M$. The root cannot fail since total wealth is at most $|E|$ and each initial wealth is one. The published market implementation extends to any fixed arity, with one more guess dimension for each extra slot before quotienting. [Original induction, market and explicit schedule](https://github.com/openai/math/blob/3344da59991ed94f9ba711217b29c1d760badfa6/preprints/Strongly-polynomial-linear-programming-November-24-2026/build/prices.tex)

| $k$ | $\delta,\gamma,\beta$ | $\eta,u,\lambda,c$ | $\epsilon,\sigma$ | $r,s$ | $g^d>r$ |
|---|---|---|---|---|---|
| 4 | $2/11,3/17,4$ | $1/16,4,1/4,2^{-14}$ | $3/1000,7/10000$ | $205001,9048$ | $d=48$ |
| 5 | $5/33,7/48,5$ | $1/16,12,1/3,2^{-18}$ | $36/10000,88/100000$ | $350001,21592$ | $d=39$ |

All strict numerical assertions can be checked over rationals. Use the fourth Taylor partial sum as a lower bound on $\ln(1+\eta)$, $-\eta-\eta^2/(2(1-\eta))$ as a lower bound on $\ln(1-\eta)$, $-c/(1-c)$ as a lower bound on $\ln(1-c)$, and $0.7\log_2(2/c)$ as an upper bound on $\ln(2/c)$ for dyadic inverse $c$. The power test is $\beta^d>r^2(1+(\beta-1)\delta)^{2d}$. [Exact parameter checks](../bench/check_price_tree.py)
