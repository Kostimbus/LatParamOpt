# Version 3.0

import argparse
from z3 import *
from sympy import primerange, isprime
from estimator import *  # provides SIS, oo, etc.
import re
import csv
import sys
from pathlib import Path
from decimal import Decimal, getcontext, ROUND_HALF_UP
import math


getcontext().prec = 80

# ------------------------- Helpers reused (trimmed) -------------------------

def resource_path(filename: str) -> str:
    if getattr(sys, "frozen", False):
        base = Path(sys._MEIPASS)  # type: ignore[attr-defined]
    else:
        base = Path(__file__).resolve().parent
    return str(base / filename)

_SEGMENTS_BY_PATH = {}

def _parse_pow2_token(tok: str) -> ArithRef:
    t = tok.strip()
    if t.startswith("2**"):
        k = int(t[3:])
        return RealVal(str(1 << k))
    return RealVal(t)

def get_log2_segments(csv_path: str):
    if csv_path in _SEGMENTS_BY_PATH:
        return _SEGMENTS_BY_PATH[csv_path]
    segs = []
    try:
        with open(csv_path, newline="") as f:
            rdr = csv.DictReader(f)
            for r in rdr:
                L = _parse_pow2_token(r["start_range"])
                U = _parse_pow2_token(r["end_range"])
                m = RealVal(r["slope"])       # accepts "num/den"
                b = RealVal(r["intercept"])   # accepts "num/den"
                segs.append((L, U, m, b))
    except FileNotFoundError:
        # Small-mode fallback: compute ≤ 0.25 error pieces on [1, 16e6]
        import math
        ln2 = math.log(2.0)
        def f(x): return math.log(x, 2.0)
        def secant(L,U):
            m = (f(U)-f(L))/(U-L); b = f(L)-m*L; return m,b
        def max_err(L,U):
            m,b = secant(L,U); xstar = max(L, min(U, 1.0/(m*ln2)))
            return f(xstar) - (m*xstar + b)
        def find_U(L,E,Umax):
            U = min(max(L+1e-9, L*1.1+1.0), Umax)
            if max_err(L,U) > E: lo,hi = L,U
            else:
                lo,hi = U, min(Umax, U*2.0)
                while hi < Umax and max_err(L,hi) <= E:
                    lo,hi = hi, min(Umax, hi*2.0)
            if hi == Umax and max_err(L,hi) <= E: return hi
            for _ in range(120):
                mid = (lo+hi)/2.0
                if max_err(L,mid) <= E: lo = mid
                else: hi = mid
                if hi-lo < 1e-13*max(1.0,lo): break
            return lo
        L, Umax, E = 1.0, 16_000_000.0, 0.25
        while L < Umax:
            U = find_U(L, E, Umax)
            if U - L < 1e-12: U = min(Umax, L + max(1e-9, L*1e-9))
            m,b = secant(L,U)
            segs.append((RealVal(str(L)), RealVal(str(U)),
                         RealVal(str(m)), RealVal(str(b))))
            L = U
    _SEGMENTS_BY_PATH[csv_path] = segs
    return segs

_LOG2_HELPERS = {}

def attach_log2_lower_bound(opt: Optimize, v: ArithRef, helper_name: str,
                            csv_path: str,
                            domain_default=(1, 16_000_000)) -> ArithRef:
    assert is_const(v) and v.num_args() == 0, "log2(...) expects a single symbol"
    vname = v.decl().name()
    key = (vname, helper_name, csv_path)
    if key in _LOG2_HELPERS:
        return _LOG2_HELPERS[key]

    y = Real(helper_name)
    _LOG2_HELPERS[key] = y

    segs = get_log2_segments(csv_path)
    vR = ToReal(v) if is_int(v) else v

    clauses = [ And(vR >= L, vR <= U, y == m*vR + b) for (L,U,m,b) in segs ]
    opt.add(Or(*clauses))

    dom = DOMAINS.get(vname, domain_default)
    opt.add(vR >= dom[0], vR <= dom[1])
    return y

def _cast_int_divisions(expr: str, z3_vars: dict) -> str:
    def is_int_name(name: str) -> bool:
        v = z3_vars.get(name)
        return v is not None and v.sort().kind() == Z3_INT_SORT

    def repl_name_num(m):
        name, num = m.group(1), m.group(2)
        if is_int_name(name):
            return f"ToReal({name}) / {num}"
        return m.group(0)

    def repl_num_name(m):
        num, name = m.group(1), m.group(2)
        if is_int_name(name):
            return f"{num} / ToReal({name})"
        return m.group(0)

    def repl_name_name(m):
        a, b = m.group(1), m.group(2)
        if is_int_name(a) and is_int_name(b):
            return f"ToReal({a}) / ToReal({b})"
        if is_int_name(a) and not is_int_name(b):
            return f"ToReal({a}) / {b}"
        if not is_int_name(a) and is_int_name(b):
            return f"{a} / ToReal({b})"
        return m.group(0)

    expr = re.sub(r"\b([A-Za-z_]\w*)\s*/\s*(\d+)\b", repl_name_num, expr)
    expr = re.sub(r"\b(\d+)\s*/\s*([A-Za-z_]\w*)\b", repl_num_name, expr)
    expr = re.sub(r"\b([A-Za-z_]\w*)\s*/\s*([A-Za-z_]\w*)\b", repl_name_name, expr)
    return expr

# parse rop exponent from any string containing "2^e"

def parse_rop_bits(text: str) -> float:
    try:
        m = re.search(r"2\^([0-9.+-]+|inf|INF|∞)", text)
        if not m:
            return 0.0
        e = m.group(1)
        if e.lower() in ("inf", "∞"):
            return float("inf")
        return float(e)
    except Exception:
        return 0.0

# Nice-print Z3 numeric values

def z3_num_to_str(val, places=4):
    if is_int_value(val):
        return str(val.as_long())
    if is_rational_value(val):
        num = Decimal(val.numerator_as_long())
        den = Decimal(val.denominator_as_long())
        dec = (num / den).quantize(Decimal(10) ** -places, rounding=ROUND_HALF_UP)
        return format(dec, 'f')
    try:
        s = val.as_decimal(places + 10)
        if s.endswith('?'):
            s = s[:-1]
        dec = Decimal(s).quantize(Decimal(10) ** -places, rounding=ROUND_HALF_UP)
        return format(dec, 'f')
    except Exception:
        return str(val)

# small primality helper (same idea as your LWE script)
SIEVE_BOUND = 2000
SMALL_PRIMES = [p for p in primerange(3, SIEVE_BOUND + 1)]
prime_vars = set()

def is_prime_sieve_expr(v: ArithRef):
    if v.sort().kind() != Z3_INT_SORT:
        raise ValueError("is_prime(var) requires an Int variable.")
    return Or(
        v == 2,
        And(
            v >= 3,
            v % 2 == 1,
            *[(v == p) | (v % p != 0) for p in SMALL_PRIMES if p > 2]
        )
    )

def is_prime_eval(var_name: str, z3_vars):
    if var_name not in z3_vars:
        raise ValueError(f"Unknown symbol '{var_name}' in is_prime({var_name})")
    prime_vars.add(var_name)
    return is_prime_sieve_expr(z3_vars[var_name])

def collect_composite_blocks(model, z3_vars):
    blocks = []
    for name in prime_vars:
        v = z3_vars[name]
        val_ast = model.eval(v, model_completion=True)
        if is_int_value(val_ast):
            val = val_ast.as_long()
        else:
            if is_rational_value(val_ast):
                num = val_ast.numerator_as_long()
                den = val_ast.denominator_as_long()
                if den != 1:
                    blocks.append((name, val_ast))
                    continue
                val = num
            else:
                continue
        if not isprime(val):
            blocks.append((name, val))
    return blocks

def estimate_bits_with_feasible_L(nv: int, qv: int) -> float:
    opt.push()
    opt.add(n == nv, q == qv)
    if opt.check() != sat:
        opt.pop()
        return float("-inf")
    mdl = opt.model()
    L_mid = mdl[L_var].as_long()
    opt.pop()

    res = SIS.estimate(SIS.Parameters(n=nv, q=qv, length_bound=L_mid, norm=norm, m=args.m))

    try:
        if "lattice" in res:
            return parse_rop_bits(res["lattice"].str())
        return parse_rop_bits(res.str())
    except Exception:
        return parse_rop_bits(res.str())

def refine_binary_search(opt, n_sym, q_sym, len_sym,
                         norm, m_opt,
                         n_lo, n_hi, q_lo, q_hi,
                         bit_sec_required, z3_vars, print_prefix=""):
    """
    Refine n, then q via binary search within the *provided* windows:
      n in [n_lo, n_hi], q in [q_lo, q_hi].
    Keeps L free (Z3 picks it). Returns (best_n, best_q).
    """
    n_var, q_var, L_var = z3_vars[n_sym], z3_vars[q_sym], z3_vars[len_sym]

    def estimate_bits_at(nv: int, qv: int) -> float:
        # Fix n,q; solve for a feasible L; then estimate
        opt.push()
        opt.add(n_var == nv, q_var == qv)
        if opt.check() != sat:
            opt.pop()
            return float("-inf")
        mdl = opt.model()
        L_mid = mdl[L_var].as_long()
        opt.pop()

        # Positional estimator call (+ optional m)
        res = SIS.estimate(SIS.Parameters(n=nv, q=qv, length_bound=L_mid, norm=norm, m=m_opt))
        
        try:
            if "lattice" in res:
                return parse_rop_bits(res["lattice"].str())
            return parse_rop_bits(res.str())
        except Exception:
            return parse_rop_bits(res.str())

    # Base refine frame: restrict the search to the *previous* windows
    opt.push()
    opt.add(n_var >= n_lo, n_var <= n_hi)
    opt.add(q_var >= q_lo, q_var <= q_hi)

    best_n, best_q = n_hi, q_hi

    # 1) Refine n in [n_lo, best_n] with q fixed at best_q
    lo, hi = n_lo, best_n
    while lo <= hi:
        mid = (lo + hi) // 2
        opt.push()
        opt.add(n_var == mid, q_var == best_q)
        ok = (opt.check() == sat)      # 'sat' here is the Z3 constant
        mdl = opt.model() if ok else None
        opt.pop()

        if sat:
            # (Optional) lazy prime blocking — if you used it elsewhere
            blocks = collect_composite_blocks(mdl, z3_vars)
            if blocks:
                for name, val in blocks:
                    opt.add(z3_vars[name] != val)
                continue  # retry same [lo, hi]

            bits = estimate_bits_at(mid, best_q)
            print(f"{print_prefix}🔍 refine n: try {mid} → ≈2^{bits:.1f}")
            if bits >= bit_sec_required:
                best_n = mid
                hi = mid - 1
            else:
                lo = mid + 1
        else:
            lo = mid + 1

    # 2) Refine q in [q_lo, best_q] with n fixed at best_n
    lo, hi = q_lo, best_q
    while lo <= hi:
        mid = (lo + hi) // 2
        opt.push()
        opt.add(n_var == best_n, q_var == mid)
        ok = (opt.check() == sat)      # 'sat' here is the Z3 constant
        mdl = opt.model() if ok else None
        opt.pop()

        if sat:
            blocks = collect_composite_blocks(mdl, z3_vars)
            if blocks:
                for name, val in blocks:
                    opt.add(z3_vars[name] != val)
                continue

            bits = estimate_bits_at(best_n, mid)
            print(f"{print_prefix}🔍 refine q: try {mid} → ≈2^{bits:.1f}")
            if bits >= bit_sec_required:
                best_q = mid
                hi = mid - 1
            else:
                lo = mid + 1
        else:
            lo = mid + 1

    opt.pop()  # drop base refine frame
    return best_n, best_q


# ------------------------- CLI -------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-bit", required=True, type=int, help="Target bit security")
parser.add_argument("-N", required=True, help="Variable name for total dimension (N = n*d)")
parser.add_argument("-d", required=True, help="Variable name for module rank d")
parser.add_argument("-n", required=True, help="Variable name for ring dimension (e.g., n)")
parser.add_argument("-q", required=True, help="Variable name for modulus (e.g., q)")
parser.add_argument("-len", dest="length", required=True, help="Variable name for length bound (e.g., L)")
parser.add_argument("-beta", required=True, type=int, choices=[0,2], help="Norm selector: 2 for euclidean, 0 for infinity (oo)")
parser.add_argument("-m", type=int, default=None, help="Length of the SIS input (m). Default: None (estimator default)")
parser.add_argument("-c", action='append', default=[], help="Z3 constraint, supports log2(...) and is_prime(x)")
parser.add_argument("-v", action='append', default=[], help="Define custom var: 'int p = q / 2' or 'real r = log2(q)'")
parser.add_argument("--enable-refine", action="store_true", help="After finding a valid solution, binary-refine n, then q within the last windows.")
parser.add_argument(
    "--bit-sec-margin",
    dest="bit_sec_margin",
    type=int,
    default=10,
    help="Security margin in bits (default: 10)"
)
parser.add_argument("-opt", "--optimize", action="append", default=[], help="Expression(s) to minimize, e.g. -opt n -opt 'log2(q)'")
parser.add_argument("--opt-priority", choices=["lex","pareto","box"], default="lex")
parser.add_argument("--enable-big", dest="enable_big", action="store_true", help="Enable big-number mode (domains up to 2**500; log2 CSV in pow2 notation)")
args = parser.parse_args()

bigMode = args.enable_big

# Base steps and ranges (aligned with your LWE script)
N_BASE_STEP = 100 if bigMode else 40
Q_BASE_STEP = 2000 if bigMode else 200

_N_BASE_STEP = N_BASE_STEP
_Q_BASE_STEP = Q_BASE_STEP

# n_range = range(41, 35001) if bigMode else range(41, 1251)
# q_range = range(300, 2**500 + 1) if bigMode else range(300, 131073)

N_MIN, N_MAX = (41, 35000) if bigMode else (41, 1250)     # inclusive
Q_MIN, Q_MAX = (300, (1 << 500) if bigMode else 131072)   # inclusive

# length domain
if bigMode:
    LENGTH_MIN, LENGTH_MAX = 50, 2**500
else:
    LENGTH_MIN, LENGTH_MAX = 50, 16_000_000

length_lb_enforced = LENGTH_MIN

LENGTH_STEP = 512 if not bigMode else 4096
length_inf_factor = 1


DEFAULT_POW2_DOMAIN = (1, 2**500) if bigMode else (1, 16_000_000)

def powers_of_two_between(lo: int, hi: int):
    lo = max(1, int(lo))
    hi = int(hi)
    vals = []
    p = 1
    # jump to first >= lo
    while p < lo:
        p <<= 1
        if p <= 0:  # overflow guard (theoretical)
            break
    while p <= hi:
        vals.append(p)
        p <<= 1
        if p <= 0:
            break
    return vals

def is_power_of_two_constraint(z3_vars, var_name, lo, hi):
    v = z3_vars[var_name]
    vals = powers_of_two_between(lo, hi)
    if not vals:
        # No power-of-two in the domain => unsat if required strictly;
        # but it's nicer to return False and let the user see UNSAT.
        return False
    if v.sort().kind() == Z3_INT_SORT:
        alts = [ v == IntVal(p) for p in vals ]
    else:
        # If user made it Real, still okay: equality to integer constants
        alts = [ v == RealVal(p) for p in vals ]
    return Or(alts)


def infer_domain_for_var(var_name, args):
    if var_name == args.n:
        return (N_MIN, N_MAX)
    if var_name == args.q:
        return (Q_MIN, Q_MAX)
    if var_name == args.len:
        return (LENGTH_MIN, LENGTH_MAX)
    # Unknown custom var: use a safe default
    return DEFAULT_POW2_DOMAIN

# For log2 helper
CSV_BASENAME = (
    "log2_big-number_mode__pow2_exp_notation__up_to_2__500_.csv" if bigMode
    else "log2_secant_segments__max_error_0_25_.csv"
)
CSV_PATH = resource_path(CSV_BASENAME)

# Known domains for helpers
# DOMAINS = {
#     args.n: (min(n_range), max(n_range)),
#     args.q: (min(q_range), max(q_range)),
#     args.length: (LENGTH_MIN, LENGTH_MAX),
# }

DOMAINS = {
    args.N: (N_MIN, N_MAX),
    args.n: (N_MIN, N_MAX),
    args.d: (1, N_MAX),
    args.q: (Q_MIN, Q_MAX),
    args.length: (LENGTH_MIN, LENGTH_MAX),
}

# Z3 vars (note: 'length' may be named 'len' by user, but we store by their chosen name)
z3_vars = {
    args.N: Int(args.N),
    args.d: Int(args.d),
    args.n: Int(args.n),
    args.q: Int(args.q),
    args.length: Int(args.length),
}

N = z3_vars[args.N]
d = z3_vars[args.d]
n = z3_vars[args.n]
q = z3_vars[args.q]
L_var = z3_vars[args.length]

# hard-coded constraints on length bound
opt = Optimize()
opt.set(priority=args.opt_priority)

# L < (q-1)/2  and  L >= 35  (cast to Real to avoid int truncation issues)
opt.add(ToReal(L_var) < (ToReal(q) - 1) / 2)
opt.add(L_var >= LENGTH_MIN)
opt.add(L_var <= LENGTH_MAX)
opt.add(N == d * n)
opt.add(N >= 41)

# Handle custom variables
if args.v:
    for definition in args.v:
        m = re.match(r"(int|real)\s+(\w+)\s*=\s*(.+)", definition)
        if m:
            var_type, var_name, expr_str = m.groups()
            var = Int(var_name) if var_type == "int" else Real(var_name)
            z3_vars[var_name] = var
            try:
                expr = eval(expr_str, {"ToReal": ToReal}, z3_vars)
                opt.add(var == expr)
            except Exception as e:
                print(f"❌ Failed to evaluate -v '{definition}': {e}")
        else:
            print(f"❌ Failed to parse -v definition: '{definition}'")

# Handle constraints (supports log2 and is_prime)
if args.c:
    for c in args.c:
        if "is_prime(" in c:
            var = re.search(r"is_prime\((\w+)\)", c).group(1)
            expr = is_prime_eval(var, z3_vars)
            opt.add(expr)
        elif "is_power_of_two(" in c:
            var = re.search(r"is_power_of_two\((\w+)\)", c).group(1)
            lo, hi = infer_domain_for_var(var, args)
            expr = is_power_of_two_constraint(z3_vars, var, lo, hi)
            if expr is False:
                print(f"⚠️ No power-of-two values for '{var}' in [{lo}, {hi}]")
                opt.add(False)  # or skip adding to just warn
            else:
                opt.add(expr)
        else:
            # bind log2 helpers then eval
            for nm in re.findall(r"log2\(\s*([A-Za-z_]\w*)\s*\)", c):
                helper = f"log{nm}"
                if helper not in z3_vars:
                    z3_vars[helper] = attach_log2_lower_bound(opt, z3_vars[nm], helper, CSV_PATH)
            c2 = re.sub(r"log2\(\s*([A-Za-z_]\w*)\s*\)", r"log\1", c)
            c3 = _cast_int_divisions(c2, z3_vars)
            try:
                expr = eval(c3, {"ToReal": ToReal}, z3_vars)
                opt.add(expr)
            except Exception as e:
                print(f"❌ Failed to parse constraint '{c}': {e}")

# Optimize objectives
if args.optimize:
    for item in args.optimize:
        # bind any log2 helpers
        for nm in re.findall(r"log2\(\s*([A-Za-z_]\w*)\s*\)", item):
            helper = f"log{nm}"
            if helper not in z3_vars:
                z3_vars[helper] = attach_log2_lower_bound(opt, z3_vars[nm], helper, CSV_PATH)
        item2 = re.sub(r"log2\(\s*([A-Za-z_]\w*)\s*\)", r"log\1", item)
        item3 = _cast_int_divisions(item2, z3_vars)
        try:
            expr = eval(item3, {"ToReal": ToReal}, z3_vars)
            opt.minimize(expr)
        except Exception as e:
            print(f"❌ Failed to parse -opt '{item}': {e}")

# ------------------------- Search -------------------------
BIT_SEC = args.bit
BIT_SEC_MARGIN = args.bit_sec_margin

N_lower_bound = N_MIN
q_lower_bound = Q_MIN

attempt = 0
acc_step_factor = 1

CSV_NOTE_PRINTED = False

no_solution = False

while True:
    attempt += 1
    opt.push()

    # current windows for n and q; L_var is free in its domain (subject to constraints)
    current_N_window_start = N_lower_bound
    current_q_window_start = q_lower_bound

    opt.add(N >= current_N_window_start, N < current_N_window_start + N_BASE_STEP)
    opt.add(q >= current_q_window_start, q < current_q_window_start + Q_BASE_STEP)

    print(f"\n🔍 Attempt {attempt}:")
    print(f"   {args.N} ∈ [{N_lower_bound}, {N_lower_bound + N_BASE_STEP})")
    print(f"   {args.q} ∈ [{q_lower_bound}, {q_lower_bound + Q_BASE_STEP})")

    if opt.check() == sat:
        model = opt.model()

        # Lazy primality blocking if needed
        blocks = collect_composite_blocks(model, z3_vars)
        if blocks:
            opt.pop()
            for name, val in blocks:
                print(f"⛔ composite candidate: {name} = {val}; blocking and retrying...")
                opt.add(z3_vars[name] != val)
            continue

        N_val = model[N].as_long()
        q_val = model[q].as_long()
        L_val = model[L_var].as_long()
        print(f"✅ Z3 solution found: {args.N}={N_val}, {args.q}={q_val}, {args.length}={L_val}")

        # Choose norm from -beta: 2 -> 2, 0 -> oo
        norm = 2 if args.beta == 2 else oo

        # Estimation

        res = SIS.estimate(SIS.Parameters(n=N_val, q=q_val, length_bound=L_val, norm=norm, m=args.m))
        

        # Try to read from keyed entry first, else the whole string
        rop_bits = 0.0
        try:
            if "lattice" in res:
                rop_bits = parse_rop_bits(res["lattice"].str())
            else:
                rop_bits = parse_rop_bits(res.str())
        except Exception:
            rop_bits = parse_rop_bits(res.str())

        print(f"Bit security: {rop_bits:.1f}")
        if math.isinf(rop_bits):
            # We want this constraint to persist across attempts,
            # so drop the per-attempt window (pop) *first*, then add the bump constraint.
            opt.pop()

            bump = LENGTH_STEP * length_inf_factor
            new_len_lb = L_val + bump
            length_lb_enforced = max(length_lb_enforced, new_len_lb + 1)

            # Tighten the global lower bound for length
            opt.add(L_var > new_len_lb)
            print(f"Infinity hardness detected (≈2^inf). Increasing {args.length} lower bound: "
                f"{args.length} > {new_len_lb} (step={LENGTH_STEP} times factor={length_inf_factor}).")

            length_inf_factor += 1
            # Try again next loop iteration with stricter length bound; keep n/q windows as-is
            continue
        else:
            # Reset the consecutive-∞ factor once we see a finite estimate
            length_inf_factor = 1


        if rop_bits >= BIT_SEC:
            # done
            if args.enable_refine and rop_bits > BIT_SEC + BIT_SEC_MARGIN:
                prev_N_start = max(current_N_window_start - N_BASE_STEP, N_MIN)
                prev_q_start = max(current_q_window_start - Q_BASE_STEP, Q_MIN)
                
                opt.pop()

                refined_N, refined_q = refine_binary_search(
                    opt, args.N, args.q, args.length, norm, args.m,
                    prev_N_start, prev_q_start,
                    N_val, q_val,
                    BIT_SEC,
                    z3_vars, print_prefix=""
                )

                # Fetch L for the refined (n,q) and re-estimate for logging
                opt.push()
                opt.add(N == refined_N, q == refined_q)
                assert opt.check() == sat
                mdl_final = opt.model()
                L_val = mdl_final[L_var].as_long()
                opt.pop()

                res2 = SIS.estimate(SIS.Parameters(n=refined_N, q=refined_q, length_bound=L_val, norm=norm, m=args.m))

                try:
                    rop_bits = parse_rop_bits(res2["lattice"].str())
                except Exception:
                    rop_bits = parse_rop_bits(res2.str())

                N_val, q_val = refined_N, refined_q
                print(f"🔧 Refined solution → {args.N}={N_val}, {args.q}={q_val}, {args.length}={L_val}")
                print(f"Bit security (after refine): ≈2^{rop_bits:.1f}")
                

            opt.push()
            opt.add(N == N_val, q == q_val)
            assert opt.check() == sat
            final_model = opt.model()

            # Print all bound values from the FINAL model (n,q forced to final)
            names_sorted = sorted(z3_vars.keys())
            for name in names_sorted:
                if name == args.N:
                    print(f"{name} = {N_val}")
                    continue
                if name == args.q:
                    print(f"{name} = {q_val}")
                    continue
                v = final_model.eval(z3_vars[name], model_completion=True)
                if v is None:
                    print(f"{name} is undefined")
                elif is_true(v) or is_false(v):
                    print(f"{name} = {v}")
                else:
                    print(f"{name} = {z3_num_to_str(v, 4)}")

            opt.pop()  # optional since we're exiting next
            break

        # else: below target → increase n and/or q
        else:
            # Hardness below target -> plan window moves using quick SAT checks
            opt.pop()  # drop the per-attempt window we just solved

            def fast_sat(N_lb, q_lb, n_shift, q_shift, n_width, q_width):
                """Check (SAT-only) if a window is feasible under all constraints."""
                opt.push()
                opt.add(n >= N_lb + n_shift, n < N_lb + n_shift + n_width)
                opt.add(q >= q_lb + q_shift, q < q_lb + q_shift + q_width)
                ok = (opt.check() == sat)
                opt.pop()
                return ok

            # 1) Try moving BOTH windows (no expansion)
            both_ok = fast_sat(
                N_lower_bound, q_lower_bound,
                N_BASE_STEP, Q_BASE_STEP,
                N_BASE_STEP, Q_BASE_STEP
            )

            if both_ok:
                #n_lower_bound += N_BASE_STEP
                N_lower_bound += N_BASE_STEP
                q_lower_bound += Q_BASE_STEP
                N_BASE_STEP = _N_BASE_STEP
                Q_BASE_STEP = _Q_BASE_STEP
                no_solution = False
                print(f"⬆️  {args.N} and {args.q}: moving both windows.")
                continue

            # 2) Else, try moving ONLY n (keep q where it is)
            only_n_ok = fast_sat(
                N_lower_bound, q_lower_bound,
                N_BASE_STEP, 0,
                N_BASE_STEP, Q_BASE_STEP
            )

            if only_n_ok:
                N_lower_bound += N_BASE_STEP
                # EXPAND the q range (interpretation: widen window size)
                N_BASE_STEP = _N_BASE_STEP
                Q_BASE_STEP *= 2
                no_solution = False
                print(f"⬆️  {args.N}: moved window; {args.q}: expanded range to {Q_BASE_STEP}.")
                continue

            # 3) Else, try moving ONLY q (keep n where it is)
            only_q_ok = fast_sat(
                N_lower_bound, q_lower_bound,
                0, Q_BASE_STEP,
                N_BASE_STEP, Q_BASE_STEP
            )

            if only_q_ok:
                q_lower_bound += Q_BASE_STEP
                # EXPAND the n range (widen window size)
                N_BASE_STEP *= 2
                Q_BASE_STEP = _Q_BASE_STEP
                no_solution = False
                print(f"⬆️  {args.q}: moved window; {args.N}: expanded range to {N_BASE_STEP}.")
                continue

            # 4) If none feasible, expand both ranges (no movement yet)
            N_BASE_STEP *= 2
            Q_BASE_STEP *= 2
            print(f"↔️  Neither move is feasible; expanding both ranges: "
                f"N_BASE_STEP={N_BASE_STEP}, Q_BASE_STEP={Q_BASE_STEP}.")
            continue

    else:
        print("❌ No solution in current range. Expanding windows.")
        opt.pop()

        no_solution = True

        if N_lower_bound > N_MAX:
            print(f"⛔ Reached {args.n} upper bound: lower_bound={N_lower_bound} > max={N_MAX}. Exiting.")
            sys.exit(1)

        if q_lower_bound > Q_MAX:
            print(f"⛔ Reached {args.q} upper bound: lower_bound={q_lower_bound} > max={Q_MAX}. Exiting.")
            sys.exit(1)

        if length_lb_enforced > LENGTH_MAX:
            print(f"⛔ Required {args.length} exceeds maximum: lower_bound={length_lb_enforced} > max={LENGTH_MAX}. Exiting.")
            sys.exit(1)

        # Otherwise, keep widening the windows as before
        acc_step_factor += 1
        N_BASE_STEP *= 2
        Q_BASE_STEP *= 2

# End
