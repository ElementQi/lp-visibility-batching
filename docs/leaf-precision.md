# Sharp precision for the balanced leaf-only schedule

A full refresh at each internal node uses logarithmic precision. If the balanced elimination reads only the normalized rows at its leaves, its general upper bound at rank proportional to dimension is quadratic in the logarithm. The second logarithm is necessary for this schedule. Exact signs of all returned leaf entries can be supplied. The obstruction applies to the rounded transcript even if upstream minor operations are exact.

## A staircase with a hidden positive direction

Let the rank be $r=2^h$ with $h\ge4$. Construct an invertible upper triangular $T_r$ and a positive vector $a_r$. Start with $T_1=[1]$, $a_1=(1)$. Write $l=r/2$, $w=\max(1,l/2)$ (at $l=1$ use $w=1$), and let $J$ be the $l$-by-$l$ matrix having ones in its first $w$ columns and zeros elsewhere. Set

$$
T_r=\begin{pmatrix} I_l&-J\\0&T_l\end{pmatrix},
\qquad a_r=(\boldsymbol 1_l,a_l/w).
$$

The first $w$ entries of $a_l$ equal one. Thus $T_ra_r=\epsilon_r e_r$, where

$$
\epsilon_r=(a_r)_r=2^{-(h-1)(h-2)/2}.
$$

Here and in the formula for the exponent we use $h\ge2$. Every coordinate of $a_r$ lies in $(0,1]$. The inverse of $T_r$ is nonnegative since $I-T_r$ is nonnegative and strictly upper triangular.

Use the two column graphs $v^\pm=-2\epsilon_r e_r\pm a_r$. For the sign-preserving, above-clipping variant below also use $v^\pm=-e_r/2\pm a_r$. Embed the flag in the original rational row matrices

$$A^\pm=[I_r\mid T_r\mid T_r v^\pm].$$

The initial echelon flag of each is precisely the displayed row order. On the first $r$ coordinates put price $2r$ and on the remainder price zero, so the effective original row matrices are $[2^{-2r}I_r\mid T_r\mid T_r v^\pm]$. We consider the first coordinate of the middle block. In the plus input the row combination $e_1^\top T_r^{-1}$ gives the middle unit row and final coordinate one. The initial block stays in $[0,1]$ because each entry of the inverse counts directed paths in an acyclic graph on $r$ vertices, at most $2^r$. In the minus input the vector $(0,-v^-,1)$ is polar and nonnegative with the same head equal to one and norm one. The true head therefore demands opposite signs. The two constructions include the same price, dimensions and initial echelon basis.

## What the balanced elimination actually reads

At every leaf, provided all previous pivots are the corresponding middle columns, the cofactor row is its original row up to common normalization. Indeed each row vanishes at the earlier middle pivots, whose prefix minor is one. That row's diagonal entry in $T_r$ is one. Other unit-magnitude coefficients come later, so the prescribed first-max tie rule selects its own middle pivot. Initial identity columns have priced magnitudes strictly less than one. All signs, including zeros, agree at every returned leaf because their only differing last-column values occur in the last row and have the same negative sign.

For center $-2\epsilon_r e_r$, every last-column leaf magnitude is at most $3\epsilon_r$. If $3\epsilon_r<2^{-E}$ then every such entry and every initial-block entry is clipped to zero. The two transcripts give the middle matrix $T_r$ and zero elsewhere. The adaptive merges on this transcript never exchange. Every dyadic interval off the nested right spine lies inside one of the identity blocks and its cross-coefficients are zero or minus one. Every dyadic interval on the right spine ends at the last middle column and has an identity graph on its own middle pivots with zero external coefficients. Hence no stabilization sees a coefficient exceeding two. This also proves the assumption on previous pivots by induction in the algorithm's depth-first order. The final chart contains only the middle identity and zero free block, independent of the true sign. Together with the algorithm's upper bound this gives $\Theta(\log^2 N)$ worst-case leaf precision for $N=2r+1$. A refresh after the final rechart reads a different original cofactor and can separate the inputs.

There is a stronger fixed-transcript version. At center $-e_r/2$ the changed entry values are $-1/2\pm\epsilon_r$, both of order one, while their absolute logarithms differ by $O(\epsilon_r)$ in base two. More precisely each log magnitude is within $8\epsilon_r$ of $-1$. When $8\epsilon_r<2^{-K}$, a source promising exact signs and log proxies of the same type as some $g$ at distance $2^{-K}$ from the true logs can give *identical full numerical replies* on this particular leaf schedule. Assign the last queried minor proxy $-1$ in both instances, keep all other queried nonzero minor logs at their common values, and round other log magnitudes to nearby rationals separately in each instance. Taking $g=f$ itself gives the type agreement for every integer affine form. The completed prefix minors are these very entries up to sign, because their remaining columns are original identity columns. The signed absolute values of the other queried nonzero minors are dyadic powers or one half. The proof allows a separate consistent extension for each input. It does not prescribe the answers of a specific implementation of the online compressor.

For this common proxy transcript the elementary absolute data are exactly those at the center, with initial identity columns clipped. The positive middle block at any merged right-spine node has terminal column $-e_{\rm last}/2$. Off the spine the terminal entries come from the direct right coupling to the last row, of magnitude at most $1/2$. Hence changing the center induces no extra basis exchanges. The pivot-free terminal sign of the true full graphs differs at the first middle pivot. This also exhibits instability under relative leaf errors of order $2^{-\Theta(\log^2 N)}$ in an entry bounded away from zero. General price-oracle complexity and algorithms that use other minor lookups are outside this lower bound.

[Exact family and physical-source regression](../bench/test_leaf_precision.py) · [Matching construction and error budget](balanced-calibration.md)
