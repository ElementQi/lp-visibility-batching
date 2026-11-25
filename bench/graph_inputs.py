"""Generate exact graph coordinates for fixed-price LP feasibility kernels."""

from fractions import Fraction as F
from pathlib import Path
import argparse
from random import Random


def mps(name):
    rows = {}
    cols = {}
    rhs = {}
    bounds = {}
    section = None
    integer = False
    marked = set()
    explicit_upper = set()
    named_sets = {}
    for line in Path(name).read_text().splitlines():
        if not line.strip() or line.startswith("*"):
            continue
        tok = line.split()
        if not line[0].isspace():
            section = tok[0]
            if section not in {"NAME", "ROWS", "COLUMNS", "RHS", "BOUNDS", "ENDATA"}:
                raise ValueError(f"unsupported MPS section {section}")
            continue
        if section == "ROWS":
            rows[tok[1]] = tok[0]
        elif section == "COLUMNS":
            if any("MARKER" in v for v in tok):
                integer = tok[-1].strip("'") == "INTORG"
                continue
            if integer:
                marked.add(tok[0])
            col = cols.setdefault(tok[0], {})
            for row, val in zip(tok[1::2], tok[2::2]):
                col[row] = col.get(row, 0) + F(val.replace("D", "E"))
        elif section == "RHS":
            if named_sets.setdefault("rhs", tok[0]) != tok[0]:
                raise ValueError("multiple RHS sets")
            for row, val in zip(tok[1::2], tok[2::2]):
                rhs[row] = F(val.replace("D", "E"))
        elif section == "BOUNDS":
            typ, set_name, key, *rest = tok
            if named_sets.setdefault("bounds", set_name) != set_name:
                raise ValueError("multiple bound sets")
            if typ in {"UP", "UI", "BV", "FX", "FR", "PL"}:
                explicit_upper.add(key)
            value = F(rest[0]) if rest else 0
            lo, hi = bounds.get(key, (0, None))
            if typ in ("UP", "UI"):
                hi = value
            elif typ in ("LO", "LI"):
                lo = value
            elif typ == "BV":
                lo, hi = 0, 1
            elif typ == "FX":
                lo = hi = value
            elif typ == "FR":
                lo = hi = None
            elif typ == "MI":
                lo = None
            elif typ == "PL":
                hi = None
            else:
                raise ValueError(typ)
            bounds[key] = lo, hi
        elif section == "RANGES":
            raise ValueError("ranges require a separate encoding")
    if marked - explicit_upper:
        raise ValueError(
            "benchmark reader requires an explicit upper bound for marked integer columns"
        )
    rows = {k: v for k, v in rows.items() if v != "N"}
    sparse = {i: {} for i in rows}
    right = {k: rhs.get(k, 0) for k in rows}
    j = 0
    for var, col in cols.items():
        lo, hi = bounds.get(var, (0, None))
        offset = lo if lo is not None else (hi if hi is not None else 0)
        for k, x in col.items():
            if k in rows:
                right[k] -= x * offset
        if lo is not None and hi == lo:
            continue
        for sign in (
            [1, -1] if lo is None and hi is None else [1 if lo is not None else -1]
        ):
            for k, x in col.items():
                if k in rows and x:
                    sparse[k][j] = sign * x
            if lo is not None and hi is not None:
                k = ("bound", j)
                sparse[k] = {j: 1}
                rows[k] = "L"
                right[k] = hi - lo
            j += 1
    for k, t in rows.items():
        if t != "E":
            assert t in ("L", "G")
            sparse[k][j] = 1 if t == "L" else -1
            j += 1
        if right[k]:
            sparse[k][-1] = -right[k]
    return [
        {(j if k == -1 else k): F(v) for k, v in row.items()} for row in sparse.values()
    ], j + 1


def svm(name, count):
    # Feasibility at a fixed objective budget adds an objective row to the soft-margin SVM.
    # Split weights and intercept into nonnegative parts; include margin slack/surplus.
    items = []
    for line in Path(name).read_text().splitlines()[:count]:
        y, *features = line.split()
        items.append(
            (int(y), {int(t.split(":")[0]) - 1: F(t.split(":")[1]) for t in features})
        )
    if count < 1 or len(items) != count or any(y not in (-1, 1) for y, _ in items):
        raise ValueError(
            "expected the requested number of binary examples labelled -1/+1"
        )
    d = max(j for _, v in items for j in v) + 1
    start = 2 * (d + 1)
    target = start + 2 * count + 1
    rows = []
    for i, (y, v) in enumerate(items):
        a = {}
        for j, x in {**v, d: F(1)}.items():
            a[2 * j] = y * x
            a[2 * j + 1] = -y * x
        a[start + i] = 1
        a[start + count + i] = -1
        a[target] = -1
        rows.append(a)
    budget = {j: 1 for j in range(2 * d)}
    budget.update({start + i: 1 for i in range(count)})
    budget[target - 1] = 1
    budget[target] = -F(count, 2)
    rows.append(budget)
    return rows, target + 1


def graph(rows, n, prices, return_order=False):
    # Original C=kernel(A). Scale each primal variable by 2^-p, so original equations acquire 2^p.
    rows = [{j: F(x) * F(2) ** prices[j] for j, x in row.items() if x} for row in rows]
    remaining = set(range(n))
    piv = []
    for t in range(len(rows)):
        candidates = [(len(rows[i]), i) for i in range(t, len(rows)) if rows[i]]
        if not candidates:
            break
        _, i = min(candidates)
        rows[t], rows[i] = rows[i], rows[t]
        row = rows[t]
        j = max(row, key=lambda j: (abs(row[j]), -j))
        assert j in remaining
        d = row[j]
        rows[t] = row = {k: x / d for k, x in row.items()}
        for i, v in enumerate(rows):
            if i == t or j not in v:
                continue
            mul = v[j]
            for k, x in row.items():
                y = v.get(k, 0) - mul * x
                if y:
                    v[k] = y
                elif k in v:
                    del v[k]
        remaining.remove(j)
        piv.append(j)
    free = sorted(remaining)
    for _ in range(10000):
        best = (F(1), -1, -1)
        for i in range(len(piv)):
            row = rows[i]
            for j, v in row.items():
                if j in remaining and abs(v) > best[0]:
                    best = (abs(v), i, j)
        _, i, j = best
        if i < 0:
            break
        row = rows[i]
        d = row[j]
        old = piv[i]
        rows[i] = row = {k: x / d for k, x in row.items()}
        for t in range(len(piv)):
            if t == i or j not in rows[t]:
                continue
            v = rows[t]
            mul = v[j]
            for k, x in row.items():
                y = v.get(k, 0) - mul * x
                if y:
                    v[k] = y
                elif k in v:
                    del v[k]
        remaining.remove(j)
        remaining.add(old)
        piv[i] = j
    else:
        raise ValueError("basis exchange budget")
    free = sorted(remaining)
    matrix = [[-rows[i].get(j, 0) for i in range(len(piv))] for j in free]
    return (matrix, tuple(free + piv)) if return_order else matrix


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("source")
    p.add_argument("tag")
    p.add_argument("destination")
    p.add_argument("--svm", type=int)
    p.add_argument("--price-seed", type=int, default=0)
    args = p.parse_args()
    rows, n = svm(args.source, args.svm) if args.svm else mps(args.source)
    rng = Random(args.price_seed)
    prices = [rng.randint(-4, 4) if args.price_seed else 0 for _ in range(n)]
    m = graph(rows, n, prices)
    assert m and m[0] and max(abs(v) for row in m for v in row) <= 1
    with open(args.destination, "x") as out:
        print(args.tag, len(m), len(m[0]), file=out)
        for row in m:
            print(*[float(x).hex() for x in row], file=out)
    print(args.tag, n, len(m), len(m[0]), flush=True)
