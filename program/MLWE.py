# Version 10.0

import argparse
from z3 import *
from sympy import primerange, isprime
from estimator import *
import re
import math
import csv
import sys
from pathlib import Path
from decimal import Decimal, getcontext, ROUND_HALF_UP


ERROR_NOISE_DISTRIBUTIONS = {
    "cbinomial": ND.CenteredBinomial,
    "dgaussian": ND.DiscreteGaussian,
    "unimod": ND.UniformMod
}

SECRET_NOISE_DISTRIBUTIONS = {
    "cbinomial": ND.CenteredBinomial,
    "binary": ND.Binary,
    "ternary": ND.Ternary,
    "sbinary": ND.SparseBinary,
    "sternary": ND.SparseTernary,
    "dgaussian": ND.DiscreteGaussian,
    "unimod": ND.UniformMod
}


getcontext().prec = 80

# LOG2
def load_log2_segments(csv_path):
    segs = []
    with open(csv_path, newline='') as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            L = parse_pow2_token(r["start_range"])
            U = parse_pow2_token(r["end_range"])
            m = RealVal(r["slope"])       # rational "num/den" is fine
            b = RealVal(r["intercept"])   # rational "num/den" is fine
            segs.append((L, U, m, b))
    return segs

_SEGMENTS_BY_PATH = {}  # cache

# LOG2 csv parser for Big Numbers
def _parse_pow2_token(tok: str) -> ArithRef:
    t = tok.strip()
    if t.startswith("2**"):
        k = int(t[3:])
        return RealVal(str(1 << k))   # exact big int -> Real
    return RealVal(t)                 # handles decimals and "num/den" rationals

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
        # Fallbacks if the CSV is missing
        if bigMode:
            # Big-mode fallback: bands by exponent span d (exact rationals)
            def band_segments(exp_max=500):
                a = 0
                while a < exp_max:
                    d = 1 if a < 128 else (2 if a < 256 else 3)
                    b_exp = min(exp_max, a + d)
                    # slope m = d / (2^a * (2^d - 1))
                    num_m = d
                    den_m = (1 << a) * ((1 << d) - 1)
                    # intercept b0 = a - d / (2^d - 1)
                    num_b = a * ((1 << d) - 1) - d
                    den_b = ((1 << d) - 1)
                    L = RealVal(str(1 << a))
                    U = RealVal(str(1 << b_exp))
                    m = RealVal(f"{num_m}/{den_m}")
                    b = RealVal(f"{num_b}/{den_b}")
                    yield (L, U, m, b)
                    a = b_exp
            segs = list(band_segments())
        else:
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

def resource_path(filename: str) -> str:
    # Works for normal runs and PyInstaller
    if getattr(sys, "frozen", False):
        base = Path(sys._MEIPASS)  # type: ignore[attr-defined]
    else:
        base = Path(__file__).resolve().parent
    return str(base / filename)

def _cast_int_divisions(expr: str, z3_vars: dict) -> str:
    """Turn Int / Int into Real division when appropriate."""

    def is_int_name(name: str) -> bool:
        v = z3_vars.get(name)
        return v is not None and v.sort().kind() == Z3_INT_SORT

    # NAME / NUMBER  -> ToReal(NAME) / NUMBER   (if NAME is Int)
    def repl_name_num(m):
        name, num = m.group(1), m.group(2)
        if is_int_name(name):
            return f"ToReal({name}) / {num}"
        return m.group(0)

    # NUMBER / NAME  -> NUMBER / ToReal(NAME)   (if NAME is Int)
    def repl_num_name(m):
        num, name = m.group(1), m.group(2)
        if is_int_name(name):
            return f"{num} / ToReal({name})"
        return m.group(0)

    # NAME1 / NAME2  -> ToReal(NAME1) / ToReal(NAME2) (if both Int)
    def repl_name_name(m):
        a, b = m.group(1), m.group(2)
        if is_int_name(a) and is_int_name(b):
            return f"ToReal({a}) / ToReal({b})"
        if is_int_name(a) and not is_int_name(b):
            return f"ToReal({a}) / {b}"
        if not is_int_name(a) and is_int_name(b):
            return f"{a} / ToReal({b})"
        return m.group(0)

    # Apply in order to avoid stepping on earlier edits
    expr = re.sub(r"\b([A-Za-z_]\w*)\s*/\s*(\d+)\b", repl_name_num, expr)
    expr = re.sub(r"\b(\d+)\s*/\s*([A-Za-z_]\w*)\b", repl_num_name, expr)
    expr = re.sub(r"\b([A-Za-z_]\w*)\s*/\s*([A-Za-z_]\w*)\b", repl_name_name, expr)
    return expr

def z3_num_to_str(val, places=4):
    # Int
    if is_int_value(val):
        return str(val.as_long())
    # Rational Real (exact num/den)
    if is_rational_value(val):
        num = Decimal(val.numerator_as_long())
        den = Decimal(val.denominator_as_long())
        dec = (num / den).quantize(Decimal(10) ** -places, rounding=ROUND_HALF_UP)
        return format(dec, 'f')
    # Algebraic/other numerics: use as_decimal then quantize
    try:
        s = val.as_decimal(places + 10)  # extra digits; may end with '?'
        if s.endswith('?'):
            s = s[:-1]
        dec = Decimal(s).quantize(Decimal(10) ** -places, rounding=ROUND_HALF_UP)
        return format(dec, 'f')
    except Exception:
        return str(val)

# --- Attack output parsing helper ---
def parse_rop_bits(text: str) -> float:
    """
    Extract the exponent e from strings like '... 2^e'.
    Handles e == 'inf' or '∞'. Returns float('inf') for infinite cost.
    Falls back to 0.0 if not found/parsable.
    """
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

def print_solution(model, z3_vars, places=4):
    for name in sorted(z3_vars.keys()):
        var = z3_vars[name]
        v = model.eval(var, model_completion=True)
        if v is None:
            print(f"{name} is undefined")
            continue
        # Bools, etc.
        if is_true(v) or is_false(v):
            print(f"{name} = {v}")
        elif v.sort().kind() in (Z3_INT_SORT, Z3_REAL_SORT):
            print(f"{name} = {z3_num_to_str(v, places)}")
        else:
            print(f"{name} = {v}")

def step_increase(factor = 2):
    return (N_BASE_STEP * factor, Q_BASE_STEP * factor, T_BASE_STEP * factor)

_SQRT_HELPERS = {}  # (var_name, helper_name) -> (r_int, sqrt_real)

def attach_sqrt_floor(opt: Optimize, v: ArithRef, helper_int: str, helper_real: str):
    """
    Create (or reuse) helpers for floor(sqrt(v)):
      r : Int  with   r*r <= v < (r+1)*(r+1)
      y : Real where  y = ToReal(r)
    Returns (r, y). Also puts y into z3_vars under helper_real if you want.
    """
    assert is_const(v) and v.num_args() == 0, "sqrt(...) expects a single symbol"
    vname = v.decl().name()
    key = (vname, helper_int)
    if key in _SQRT_HELPERS:
        return _SQRT_HELPERS[key]

    r = Int(helper_int)     # floor(sqrt(v))
    y = ToReal(r)           # Real view (use this in linear Real arithmetic)

    # Domain: v must be non-negative for sqrt
    if is_int(v):
        opt.add(v >= 0)
    else:
        # if Real and you truly allow it, still require v >= 0
        opt.add(v >= 0)

    # Pin r to the correct integer interval: r^2 <= v < (r+1)^2
    # (This is NIA, same as your previous approach.)
    opt.add(r*r <= v)
    opt.add((r + 1)*(r + 1) > v)

    _SQRT_HELPERS[key] = (r, y)
    return r, y

def handle_sqrt_in_str(c_str: str, z3_vars: dict, opt: Optimize) -> str:
    """
    - Creates helpers for each sqrt(NAME)
       * Int helper: __sqrt_floor_NAME
       * Real alias: sqrtNAME := ToReal(__sqrt_floor_NAME)
    - Rewrites 'sqrt(NAME)' -> 'sqrtNAME'
    - Adds ToReal casts around * and / near sqrt terms
    Returns the rewritten string.
    """
    names = re.findall(r"sqrt\(\s*([A-Za-z_]\w*)\s*\)", c_str)
    for nm in names:
        if nm not in z3_vars:
            raise ValueError(f"Unknown symbol '{nm}' in sqrt({nm})")
        int_h = f"__sqrt_floor_{nm}"
        real_h = f"sqrt{nm}"
        if real_h not in z3_vars:
            r_int, y_real = attach_sqrt_floor(opt, z3_vars[nm], int_h, real_h)
            # store both if you want; the Real alias is what expressions will use
            z3_vars[int_h] = r_int
            z3_vars[real_h] = y_real

    # Replace sqrt(name) with the Real alias
    c2 = re.sub(r"sqrt\(\s*([A-Za-z_]\w*)\s*\)", r"sqrt\1", c_str)

    c2 = _cast_int_divisions(c2, z3_vars)

    return c2

# Avoid duplicating helpers
_LOG2_HELPERS = {}  # (var_name, helper_name, csv_path) -> RealRef


def refine_binary_search(opt, N_var, q_var, t_var, N_val, q_val, t_val,
                         last_N, last_q, last_t,
                         N_overkill, q_overkill, t_overkill,
                         bit_sec_required, margin,
                         z3_vars, attack_order, print_prefix=""):

    refined_result = (N_val, q_val, t_val)

    def secure_enough(sec):
        return all(sec.get(atk, 0) >= bit_sec_required for atk in ["usvp", "dual", "dual_hybrid", "bdd", "bkw"]) and \
               sec.get("bkw", 0) >= bit_sec_required and \
               sec.get("arora-gb", 0) >= bit_sec_required

    def too_secure(sec):
        return any(sec.get(atk, 0) > bit_sec_required + margin for atk in ["usvp", "dual", "dual_hybrid", "bdd", "bkw"]) or \
               sec.get("bkw", 0) > bit_sec_required + margin or \
               sec.get("arora-gb", 0) > bit_sec_required + margin

    def lwe_estimate(n_try, q_try, t_try):
        lwe_params = LWE.Parameters(N=n_try, q=q_try, Xs=Xs_inst, Xe=Xe_dist(t_try))
        return LWE.estimate(lwe_params)

    param_bounds = {
        N_var: (last_N, N_val) if N_overkill else (N_val, N_val+1),
        q_var: (last_q, q_val) if q_overkill else (q_val, q_val+1),
        t_var: (last_t, t_val) if t_overkill else (t_val, t_val+1)
    }

    for var, (lo, hi) in param_bounds.items():
        if hi - lo <= 1:
            continue

        low, high = lo, hi
        best_found = None

        while high - low > 1:
            mid = (low + high) // 2

            opt.push()
            # opt.add(z3_vars[var] == mid)
            opt.add(z3_vars[var] <= mid, z3_vars[var] >= low)
            opt.add(z3_vars[N_var] == N_val if var != N_var else True)
            opt.add(z3_vars[q_var] == q_val if var != q_var else True)
            opt.add(z3_vars[t_var] == t_val if var != t_var else True)


            if opt.check() == sat:
                print(f"{print_prefix}🔍 Binary search testing {var}={mid}")
                model = opt.model()

                blocks = collect_composite_blocks(model)
                if blocks:
                    opt.pop()
                    for name, val in blocks:
                        print(f"{print_prefix}   ⛔ composite during refine: {name} = {val}; blocking and retrying...")
                        opt.add(z3_vars[name] != val)
                    continue  # retry with same [low, high)

                est_n = model[z3_vars[N_var]].as_long()
                est_q = model[z3_vars[q_var]].as_long()
                t_ast = model[z3_vars[t_var]]
                if is_int_value(t_ast):
                    est_t = t_ast.as_long()
                elif is_rational_value(t_ast):
                    est_t = t_ast.numerator_as_long() / t_ast.denominator_as_long()
                else:
                    s = t_ast.as_decimal(20)
                    if s.endswith('?'):
                        s = s[:-1]
                    est_t = float(s)

                result = lwe_estimate(est_n, est_q, est_t)
                security = {}
                for atk in attack_order:
                    try:
                        #rop = float(re.search(r'rop: ≈2\^([0-9.]+)', result[atk].str()).group(1))
                        rop = parse_rop_bits(result[atk].str())
                        security[atk] = rop
                    except:
                        security[atk] = 0

                print(f"{print_prefix}   Tested ({est_n}, {est_q}, {est_t}) -> {security}")

                if secure_enough(security):
                    best_found = (est_n, est_q, est_t)
                    high = mid
                else:
                    low = mid
            else:
                print(f"{print_prefix}   No SAT for {var}={mid}, increasing lower bound")
                low = mid
            opt.pop()

        if best_found:
            refined_result = best_found

    return refined_result

def handle_log2_in_str(c_str: str, z3_vars: dict, opt: Optimize, csv_path: str) -> str:
    """Bind log2(VAR) helpers and rewrite text to use them."""
    # 1) create helpers for every log2(NAME)
    names = re.findall(r"log2\(\s*([A-Za-z_]\w*)\s*\)", c_str)
    for nm in names:
        if nm not in z3_vars:
            raise ValueError(f"Unknown symbol '{nm}' in log2({nm})")
        helper = f"log{nm}"
        if helper not in z3_vars:
            z3_vars[helper] = attach_log2_lower_bound(opt, z3_vars[nm], helper, CSV_PATH)

    # 2) rewrite log2(name) -> logname
    c2 = re.sub(r"log2\(\s*([A-Za-z_]\w*)\s*\)", r"log\1", c_str)

    c2 = _cast_int_divisions(c2, z3_vars)

    return c2


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

# ------------------------- CLI ARGUMENT PARSING -------------------------
parser = argparse.ArgumentParser()
#parser.add_argument("name", required=True, help="Problem name, e.g. MLWE")
parser.add_argument("-bit", required=True, type=int, help="Target bit security")
parser.add_argument("-N", required=True, help="Variable name for total MLWE dimension N")
parser.add_argument("-n", required=True, help="Variable name for ring dimension n")
parser.add_argument("-d", required=True, help="Variable name for module rank d") 
parser.add_argument("-q", required=True, help="Variable name for modulus (e.g., q)")
parser.add_argument("-t", required=True, help="Variable name for noise (e.g., t)")
parser.add_argument("--enable-refine", action="store_true", help="Enable the refinement of the parameters")
parser.add_argument("-c", action='append', default=[], help="Z3 constraint (e.g. 'N % 2 == 0', 'is_prime(q)', 'log2(N) < 10')")
parser.add_argument("-v", action='append', default=[], help="Z3 constraint (e.g. ")
parser.add_argument("-opt", "--optimize", action="append", default=[],
    help="Variable or linear expression to minimize, e.g. -opt N, -opt 'q/10', -opt 'log2(v)'")
parser.add_argument("--opt-priority", choices=["lex","pareto","box"], default="lex")
parser.add_argument(
    "--enable-big",
    dest="enable_big",
    action="store_true",
    help="Enable big-number mode for parameters and log2: CSV ranges as 2**k, domain up to 2**500."
)
parser.add_argument(
    "--error-distribution", 
    dest="error_distribution", 
    default="cbinomial", 
    choices=list(ERROR_NOISE_DISTRIBUTIONS.keys()),
    help="Error distribution (Xe): cbinomial | dgaussian | unimod"
)
parser.add_argument(
    "--secret-distribution", 
    dest="secret_distribution", 
    default="cbinomial", 
    choices=list(SECRET_NOISE_DISTRIBUTIONS.keys()),
    help="Secret distribution (Xs), key from SECRET_NOISE_DISTRIBUTIONS: cbinomial | binary | sbinary | sternary"
)
parser.add_argument(
    "-s", 
    type=int,
    default=None, 
    help="Optional parameter for the secret distribution"
)
parser.add_argument(
    "--bit-sec-margin",
    dest="bit_sec_margin",
    type=int,
    default=10,
    help="Security margin in bits (default: 10)"
)
args = parser.parse_args()
bigMode = args.enable_big

DEFAULT_POW2_DOMAIN = (1, 2**500) if bigMode else (1, 16_000_000)

def infer_domain_for_var(var_name, args):
    if var_name == args.N:
        return (N_MIN, N_MAX)
    if var_name == args.n:
        # ring dimension n: at least 40; upper-bounded indirectly by N via N = d*n
        return (41, N_MAX)
    if var_name == args.d:
        # module rank d: unbounded a priori; cap by max N for safety
        return (1, N_MAX)
    if var_name == args.q:
        return (Q_MIN, Q_MAX)
    if var_name == args.t:
        return (T_MIN, T_MAX)
    # Unknown custom var: use a safe default
    return DEFAULT_POW2_DOMAIN


CSV_BASENAME = (
    "log2_big-number_mode__pow2_exp_notation__up_to_2__500_.csv"
    if bigMode
    else "log2_secant_segments__max_error_0_25_.csv"
)
CSV_PATH = resource_path(CSV_BASENAME)


# ------------------------- SETUP -------------------------
BIT_SEC = args.bit
BIT_SEC_MARGIN = args.bit_sec_margin

sec_key = args.secret_distribution.strip().lower()
sec_param = args.s

if sec_key in ("binary", "ternary") :
    if sec_param is not None:
        print(f"Secret distribution '{sec_key}' ignores -s (got {sec_param}); proceeding without it.")
    Xs_inst = SECRET_NOISE_DISTRIBUTIONS[sec_key]
elif sec_key in ("sbinary", "sternary"):
    if sec_param is None:
        print(f"Secret distribution '{sec_key}' got {sec_param} for the distribution parameter; proceeding with 16.")
        sec_param = 16
    elif sec_param < 4:
        print(f"{sec_key}: provided -s={sec_param} < 4; using 4 instead.")
        sec_param = 4
    Xs_inst = SECRET_NOISE_DISTRIBUTIONS[sec_key](sec_param)
elif sec_key == "cbinomial":
    if sec_param is None:
        print(f"Secret distribution '{sec_key}' got {sec_param} for the distribution parameter; proceeding with 2.")
        sec_param = 2
    elif sec_param < 1:
        print(f" cbinomial: provided -s={sec_param} < 1; using 1 instead.")
        sec_param = 1
    Xs_inst = SECRET_NOISE_DISTRIBUTIONS[sec_key](sec_param)
elif sec_key == "dgaussian":
    if sec_param is None:
        print(f"Secret distribution '{sec_key}' got {sec_param} for the distribution parameter; proceeding with 1.5.")
        sec_param = 1.5
    elif sec_param < 0.5:
        print(f" dgaussian: provided -s={sec_param} < 0.5; using 0.5 instead.")
        sec_param = 0.5
    Xs_inst = SECRET_NOISE_DISTRIBUTIONS[sec_key](sec_param)
elif sec_key == "unimod":
    if sec_param is None:
        print(f"Secret distribution '{sec_key}' got {sec_param} for the distribution parameter; proceeding with 5.")
        sec_param = 5
    elif sec_param < 2:
        print(f" unimod: provided -s={sec_param} < 2; using 2 instead.")
        sec_param = 2
    Xs_inst = SECRET_NOISE_DISTRIBUTIONS[sec_key](sec_param)
else:
    raise SystemExit(f"Unknown secret distribution '{sec_key}'.")


err_key = args.error_distribution.strip().lower()
Xe_dist = ERROR_NOISE_DISTRIBUTIONS[err_key]
# Dynamically create Z3 Int variables
t_is_real = (err_key == "dgaussian")

N_BASE_STEP = (
    100
    if bigMode
    else 30
)

Q_BASE_STEP = (
    2000
    if bigMode
    else 200
)

if t_is_real:
    T_BASE_STEP = (
        0.5
        if bigMode
        else 0.5
    )
else:
    T_BASE_STEP = (
        1
        if bigMode
        else 1
    )

_N_BASE_STEP, _Q_BASE_STEP, _T_BASE_STEP = N_BASE_STEP, Q_BASE_STEP, T_BASE_STEP
var_names = {"N": args.N, "n": args.n, "d": args.d, "q": args.q, "t": args.t}


z3_vars = {
    args.N: Int(args.N),
    args.n: Int(args.n),
    args.d: Int(args.d),
    args.q: Int(args.q),
    args.t: Real(args.t) if t_is_real else Int(args.t),
}

N = z3_vars[args.N]
d = z3_vars[args.d]
n_ring = z3_vars[args.n]
q = z3_vars[args.q]
t = z3_vars[args.t]


N_MIN, N_MAX = (41, 35000) if bigMode else (41, 1250)     # inclusive
Q_MIN, Q_MAX = (300, (1 << 500) if bigMode else 131072)   # inclusive
if t_is_real:
    #t_range = range(0.5, 8.0) if bigMode else range(0.5, 3.5)
    T_MIN, T_MAX = (0.5, 8.0) if bigMode else (0.5, 3.5)
else:
    if err_key == "cbinomial" or err_key == "unimod":
        T_MIN, T_MAX = (1, 33) if bigMode else (1, 11)


DOMAINS = {
    args.N: (N_MIN, N_MAX),
    args.q: (Q_MIN, Q_MAX),
    args.t: (T_MIN, T_MAX),
}

LOG2_DOMAIN_DEFAULT = (1, 1 << 500) if bigMode else (1, 16_000_000)

def attach_log2_lower_bound(opt: Optimize, v: ArithRef, helper_name: str,
                            csv_path: str,
                            domain_default=LOG2_DOMAIN_DEFAULT) -> ArithRef:
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


# Initialize Z3 optimizer
opt = Optimize()


# Handle the refine parameter
refine = args.enable_refine

# Prime list for Normal Mode
#prime_list = list(primerange(300, 131073))


# Add CLI custom variables
if args.v:
    for definition in args.v:
        # Parse variable declaration
        match = re.match(r"(int|real)\s+(\w+)\s*=\s*(.+)", definition)
        if match:
            var_type, var_name, expr_str = match.groups()
            # Define variable type
            if var_type == "int":
                var = Int(var_name)
            elif var_type == "real":
                var = Real(var_name)
            else:
                print(f"❌ Unknown type '{var_type}'")
                continue
            # Add to z3_vars for later reference
            z3_vars[var_name] = var
            try:
                expr = eval(expr_str, {}, z3_vars)
                opt.add(var == expr)
            except Exception as e:
                print(f"❌ Failed to evaluate expression '{expr_str}': {e}")
        else:
            print(f"❌ Failed to parse variable definition: '{definition}'")


def handle_constraint(c_str, z3_vars, opt):
    # 1) If it contains log2(VAR), attach helper(s) and rewrite
    def bind(name):
        if name not in z3_vars:
            raise ValueError(f"Unknown symbol '{name}' in log2({name})")
        helper = f"log{name}"
        if helper not in z3_vars:
            z3_vars[helper] = attach_log2_lower_bound(opt, z3_vars[name], helper, CSV_PATH)
        return helper

    for nm in re.findall(r"log2\(\s*([A-Za-z_]\w*)\s*\)", c_str):
        bind(nm)

    # Rewrite "log2(v)" -> "logv"
    c2 = re.sub(r"log2\(\s*([A-Za-z_]\w*)\s*\)", r"log\1", c_str)

    c2 = _cast_int_divisions(c2, z3_vars)

    # 3) Evaluate with a minimal safe env (no Python builtins)
    env = {"ToReal": ToReal}
    env.update(z3_vars)

    expr = eval(c2, {}, env)
    opt.add(expr)


SIEVE_BOUND = 1000
SMALL_PRIMES = [p for p in primerange(3, SIEVE_BOUND + 1)]
prime_vars = set()

def is_prime_sieve_expr(v: ArithRef):
    if v.sort().kind() != Z3_INT_SORT:
        raise ValueError("is_prime(var) requires an Int variable.")
    # Or(v == 2, And(v >= 3, v % 2 == 1, ∧ v % p != 0))
    #return Or(v == 2, And(v >= 3, v % 2 == 1, *[v % p != 0 for p in SMALL_PRIMES]))
    return Or(
        v == 2,
        And(
            v >= 3,
            v % 2 == 1,
            *[(v == p) | (v % p != 0) for p in SMALL_PRIMES if p > 2]
        )
    )

def is_prime_eval(var_name: str):
    if var_name not in z3_vars:
        raise ValueError(f"Unknown symbol '{var_name}' in is_prime({var_name})")
    prime_vars.add(var_name)
    return is_prime_sieve_expr(z3_vars[var_name])

def collect_composite_blocks(model):
    """Return list of (var_name, value) that violate primality under the model."""
    blocks = []
    for name in prime_vars:
        v = z3_vars[name]
        val_ast = model.eval(v, model_completion=True)
        if is_int_value(val_ast):
            val = val_ast.as_long()
        else:
            # Try to coerce rationals safely
            if is_rational_value(val_ast):
                num = val_ast.numerator_as_long()
                den = val_ast.denominator_as_long()
                if den != 1:
                    # Non-integer -> cannot be prime
                    blocks.append((name, val_ast))
                    continue
                val = num
            else:
                # algebraic / non-ground: skip
                continue
        # Lazy primality check
        if not isprime(val):
            blocks.append((name, val))
    return blocks


# Add CLI constraints (with log2 handling)
if args.c:
    for c in args.c:
        if "is_prime(" in c:
            # var = re.search(r"is_prime\((\w+)\)", c).group(1)
            # expr = Or([z3_vars[var] == p for p in prime_list])
            # opt.add(expr)
            var = re.search(r"is_prime\((\w+)\)", c).group(1)
            expr = is_prime_eval(var)
            #expr = Or([z3_vars[var] == p for p in prime_list])
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
            try:
                c2 = handle_log2_in_str(c, z3_vars, opt, CSV_PATH)
                c3 = handle_sqrt_in_str(c2, z3_vars, opt)
                expr = eval(c3, {"ToReal": ToReal}, z3_vars)
                #expr = eval(c, {}, z3_vars)
                opt.add(expr)
            except Exception as e:
                print(f"❌ Failed to parse constraint '{c}': {e}")


opt.add(N == d * n_ring)
opt.add(N >= 40)

opt.set(priority=args.opt_priority)
goals = []


if args.optimize:
    for item in args.optimize:
        try:
            # same log2 handling as constraints
            item2 = handle_log2_in_str(item, z3_vars, opt, CSV_PATH)
            item3 = handle_sqrt_in_str(item2, z3_vars, opt)
            expr = eval(item3, {"ToReal": ToReal}, z3_vars)
            goals.append(opt.minimize(expr))
        except Exception as e:
            print(f"❌ Failed to parse -opt '{item}': {e}")


# ------------------------- PARAMETER SEARCH -------------------------
attack_order = ["usvp", "dual", "dual_hybrid", "bkw", "arora-gb", "bdd"]
N_lower_bound = N_MAX
q_lower_bound = Q_MIN
if err_key == "dgaussian":
    t_lower_bound = 0.5
elif err_key == "unimod":
    t_lower_bound = 2
elif err_key == "cbinomial":
    t_lower_bound = 1

attempt = 0
acc_step_factor = 1
last_N = N_lower_bound
last_q = q_lower_bound
last_t = t_lower_bound

no_solution = False


while True:
    attempt += 1
    opt.push()

    # Define current parameter ranges
    opt.add(N >= N_lower_bound, N < N_lower_bound + N_BASE_STEP)
    opt.add(q >= q_lower_bound, q < q_lower_bound + Q_BASE_STEP)
    opt.add(t >= t_lower_bound, t < t_lower_bound + T_BASE_STEP)

    print(f"\n🔍 Attempt {attempt}:")
    print(f"   {args.N} ∈ [{N_lower_bound}, {N_lower_bound + N_BASE_STEP})")
    print(f"   {args.q} ∈ [{q_lower_bound}, {q_lower_bound + Q_BASE_STEP})")
    print(f"   {args.t} ∈ [{t_lower_bound}, {t_lower_bound + T_BASE_STEP})")


    if opt.check() == sat:
        if(acc_step_factor > 1):
            acc_step_factor = 1
        model = opt.model()

        # Lazy primality check (blocking) in main loop
        blocks = collect_composite_blocks(model)
        if blocks:
            opt.pop()
            for name, val in blocks:
                print(f"⛔ composite candidate: {name} = {val}; blocking and retrying...")
                opt.add(z3_vars[name] != val)
            no_solution = False
            continue  # retry with same bounds

        N_val = model[N].as_long()
        q_val = model[q].as_long()
        #t_val = model[t].as_long()
        t_ast = model[t]
        if is_int_value(t_ast):
            t_val = t_ast.as_long()
        elif is_rational_value(t_ast):
            t_val = t_ast.numerator_as_long() / t_ast.denominator_as_long()
        else:
            s = t_ast.as_decimal(20)
            if s.endswith('?'):
                s = s[:-1]
            t_val = float(s)
        print(f"✅ Z3 solution found: {args.N}={N_val}, {args.q}={q_val}, {args.t}={t_val}")

        lwe_params = LWE.Parameters(n=N_val, q=q_val, Xs=Xs_inst, Xe=Xe_dist(t_val))
        result = LWE.estimate(lwe_params)

        security = {}
        for atk in attack_order:
            if atk in result:
                try:
                    #rop = float(re.search(r'rop: ≈2\^([0-9.]+)', result[atk].str()).group(1))
                    rop = parse_rop_bits(result[atk].str())
                    security[atk] = rop
                except Exception:
                    security[atk] = 0
            else:
                security[atk] = 0

        print(f"Bit security: {security}")

        increase_N = any(security[atk] < BIT_SEC for atk in ["usvp", "dual", "dual_hybrid", "bdd", "bkw"])
        increase_q_candidate = security.get("bkw", 0) < BIT_SEC
        increase_t = security.get("arora-gb", 0) < BIT_SEC or (any(security[atk] < BIT_SEC for atk in ["usvp", "dual", "dual_hybrid", "bdd", "bkw"]) and (err_key == "unimod" and t_val < 5))

        overkill_N = any(security[atk] > BIT_SEC + BIT_SEC_MARGIN for atk in ["usvp", "dual", "dual_hybrid", "bdd", "bkw"])
        overkill_q = security.get("bkw", 0) > BIT_SEC + BIT_SEC_MARGIN
        overkill_t = security.get("arora-gb", 0) > BIT_SEC + BIT_SEC_MARGIN

        if not (increase_N or increase_q_candidate or increase_t):
            print("🎉 All attacks passed threshold!")
            if(refine and (overkill_N or overkill_q or overkill_t)):
                opt.pop()
                refined_N, refined_q, refined_t = refine_binary_search(
                    opt, args.N, args.q, args.t,
                    N_val, q_val, t_val,
                    last_N, last_q, last_t,
                    overkill_N, overkill_q, overkill_t,
                    BIT_SEC, BIT_SEC_MARGIN,
                    z3_vars, attack_order)
                print_solution(model, z3_vars)
                print(f"✅ Refined params: {args.N}={refined_N}, {args.q}={refined_q}, {args.t}={refined_t}")
            elif (not refine):
                print_solution(model, z3_vars)
            break  # Stop if you want, or comment this out to keep searching tighter

        opt.pop()
        if increase_N:
            if(no_solution):
                N_lower_bound = N_val
                N_BASE_STEP = _N_BASE_STEP
            else:     
                N_lower_bound += N_BASE_STEP
            last_N = N_val
            print(f"⬆️  Increasing {args.N}")
        if increase_q_candidate:
            if q_val <= 3000:
                if(no_solution):
                    q_lower_bound = q_val
                    Q_BASE_STEP = _Q_BASE_STEP
                else:
                    q_lower_bound += Q_BASE_STEP
                last_q = q_val
                print(f"⬆️  Increasing {args.q}")
            else:
                # Don’t shift the lower bound; just widen the window so other vars/constraints can use larger q
                old_q_end = q_lower_bound + Q_BASE_STEP
                Q_BASE_STEP += Q_BASE_STEP
                new_q_end = q_lower_bound + Q_BASE_STEP
                print(f"🪄 BKW < target but {args.q}={q_val} > 3000 → keeping lower bound, widening {args.q} range:")
                print(f"   {args.q} ∈ [{q_lower_bound}, {old_q_end}) → [{q_lower_bound}, {new_q_end})")
        if increase_t:
            if(no_solution):
                t_lower_bound = t_val
                T_BASE_STEP = _T_BASE_STEP
            else:
                t_lower_bound += T_BASE_STEP
            last_t = t_val
            print(f"⬆️  Increasing {args.t}")
        if (increase_N or increase_q_candidate or increase_t):
            no_solution = False

    else:
        print("❌ No solution in current range. Expanding.")
        opt.pop()
        no_solution = True
        acc_step_factor += 1

        if N_lower_bound > N_MAX:
            print(f"⛔ Reached {args.n} upper bound: lower_bound={N_lower_bound} > max={N_MAX}. Exiting.")
            sys.exit(1)

        if q_lower_bound > Q_MAX:
            print(f"⛔ Reached {args.q} upper bound: lower_bound={q_lower_bound} > max={Q_MAX}. Exiting.")
            sys.exit(1)

        if t_lower_bound > T_MAX:
            print(f"⛔ Required {args.t} exceeds maximum: lower_bound={t_lower_bound} > max={T_MAX}. Exiting.")
            sys.exit(1)
        N_BASE_STEP, Q_BASE_STEP, T_BASE_STEP = step_increase(acc_step_factor)









