"""Digital Multi-Load Beam Structural Analysis and Simulation System.

Run with Python 3.11+; Tkinter is the only GUI dependency. Importing this
module does not open a window, so the calculation functions can be tested.
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
import math
import tkinter as tk
from tkinter import messagebox, ttk


APP_TITLE = "Digital Multi-Load Beam Structural Analysis and Simulation System"
NAME = "Emman James P. Faderon"
APP_DESCRIPTION = (
    "A Python-based program that analyzes a structural beam under different loading conditions, "
    "such as point loads, uniformly distributed loads, triangular loads, and applied moments. "
    "The program allows users to select different types of supports and automatically calculates "
    "the support reactions, shear force, bending moment, slope, and deflection. It also provides "
    "a digital visualization of the beam and its deformed shape to help users understand how "
    "loads affect a structural beam."
)


@dataclass(frozen=True)
class Support:
    kind: str
    x: float


@dataclass(frozen=True)
class Load:
    kind: str
    a: float
    value: float
    b: float = 0.0
    end_value: float = 0.0


@dataclass(frozen=True)
class Reaction:
    support: Support
    force: float
    moment: float


def finite_number(value: object, label: str) -> float:
    """Parse numbers once at the input boundary and reject NaN/infinity."""
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be a number.") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite.")
    return result


def stiffness_value(mode: str, direct: object, modulus: object, inertia: object) -> float:
    """Read only the active constant stiffness inputs; E uses GPa and I uses mm⁴."""
    if mode == "Constant EI":
        ei = finite_number(direct, "EI")
    elif mode == "E and I":
        e = finite_number(modulus, "E")
        i = finite_number(inertia, "I")
        if e <= 0 or i <= 0:
            raise ValueError("E and I must each be greater than zero.")
        # Separate exponents avoid unnecessary intermediate product overflow.
        em, ee = math.frexp(e)
        im, ie = math.frexp(i)
        fm, fe = math.frexp(1e-6)
        try:
            ei = math.ldexp(em * im * fm, ee + ie + fe)
        except OverflowError as exc:
            raise ValueError("E × I exceeds floating-point range.") from exc
    else:
        raise ValueError("Choose Constant EI or E and I.")
    if not math.isfinite(ei) or ei <= 0:
        raise ValueError("EI must be finite and greater than zero; check the active stiffness inputs.")
    return ei


def validate_model(length: float, ei: float, supports: list[Support],
                   loads: list[Load]) -> None:
    for value, label in ((length, "Length"), (ei, "EI")):
        if finite_number(value, label) <= 0:
            raise ValueError(f"{label} must be greater than zero.")
    positions: set[float] = set()
    for support in supports:
        if support.kind not in ("Pinned", "Roller", "Fixed"):
            raise ValueError("Choose Pinned, Roller, or Fixed support.")
        x = finite_number(support.x, "Support position")
        if not 0 <= x <= length:
            raise ValueError("Support positions must be within the beam.")
        if x in positions:
            raise ValueError("A support already exists there; edit it instead.")
        positions.add(x)
        if support.kind == "Fixed" and x not in (0, length):
            raise ValueError("Fixed supports must be at a beam end.")
    if not any(s.kind == "Fixed" for s in supports) and len(positions) < 2:
        raise ValueError("Unstable beam: use a fixed end or two distinct supports.")
    for load in loads:
        if load.kind not in ("Point", "Moment", "UDL", "Linear"):
            raise ValueError("Choose Point, Moment, UDL, or Linear load.")
        a = finite_number(load.a, "Load position")
        finite_number(load.value, "Load value")
        if not 0 <= a <= length:
            raise ValueError("Load positions must be within the beam.")
        if load.kind in ("UDL", "Linear"):
            b = finite_number(load.b, "Load end")
            finite_number(load.end_value, "End intensity")
            if not a < b <= length:
                raise ValueError("Distributed loads need 0 ≤ start < end ≤ length.")


def distributed(load: Load) -> bool:
    return load.kind in ("UDL", "Linear")


def intensity(load: Load, x: float) -> float:
    if load.kind == "UDL":
        return load.value
    return load.value + (load.end_value - load.value) * ((x - load.a) / (load.b - load.a))


def load_resultant(load: Load) -> tuple[float, float]:
    """Return downward force and its first moment about x=0."""
    if load.kind == "Point":
        return load.value, load.value * load.a
    if load.kind == "Moment":
        return 0.0, 0.0
    h = load.b - load.a
    q0, q1 = load.value, intensity(load, load.b)
    force = h * (q0 + q1) / 2
    return force, load.a * force + h * h * (q0 + 2 * q1) / 6


def element_stiffness(h: float, ei: float) -> list[list[float]]:
    factor = ei / h / h / h
    return [[factor * v for v in row] for row in (
        (12, 6*h, -12, 6*h),
        (6*h, 4*h*h, -6*h, 2*h*h),
        (-12, -6*h, 12, -6*h),
        (6*h, 2*h*h, -6*h, 4*h*h))]


def consistent_load(h: float, q0: float, q1: float) -> list[float]:
    return [-h*(7*q0+3*q1)/20, -h*h*(3*q0+2*q1)/60,
            -h*(3*q0+7*q1)/20, h*h*(2*q0+3*q1)/60]


def dense_solve(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Solve a symmetric positive system with diagonal scaling and pivoting."""
    n = len(rhs)
    if not n:
        return []
    if any(matrix[i][i] <= 0 or not math.isfinite(matrix[i][i]) for i in range(n)):
        raise ValueError("The model is unstable or its numerical scale is too extreme.")
    scale = [1 / math.sqrt(matrix[i][i]) for i in range(n)]
    original = [[matrix[i][j] * scale[i] * scale[j] for j in range(n)] for i in range(n)]
    target = [rhs[i] * scale[i] for i in range(n)]
    rows = [original[i][:] + [target[i]] for i in range(n)]
    if any(not math.isfinite(v) for row in rows for v in row):
        raise ValueError("The model exceeds floating-point range; use consistent scales.")
    for col in range(n):
        pivot = max(range(col, n), key=lambda i: abs(rows[i][col]))
        if abs(rows[pivot][col]) < 1e-12:
            raise ValueError("The model is unstable or too ill-conditioned. Check nearby nodes and scales.")
        rows[col], rows[pivot] = rows[pivot], rows[col]
        for row in range(col + 1, n):
            factor = rows[row][col] / rows[col][col]
            rows[row][col] = 0.0
            for j in range(col + 1, n + 1):
                rows[row][j] -= factor * rows[col][j]
    answer = [0.0] * n
    for i in range(n - 1, -1, -1):
        answer[i] = (rows[i][n] - math.fsum(rows[i][j] * answer[j]
                     for j in range(i + 1, n))) / rows[i][i]
    if any(not math.isfinite(v) for v in answer):
        raise ValueError("The calculation overflowed; check model scales.")
    for i in range(n):
        terms = [original[i][j] * answer[j] for j in range(n)]
        error = abs(math.fsum(terms) - target[i])
        bound = math.fsum(abs(v) for v in terms) + abs(target[i])
        if error > 1e-8 * max(bound, 1e-300):
            raise ValueError("The numerical solution did not pass its residual check.")
    return [answer[i] * scale[i] for i in range(n)]


def polynomial_value(coefficients: list[float], x: float) -> float:
    result = 0.0
    for coefficient in reversed(coefficients):
        result = result * x + coefficient
    return result


def polynomial_roots(coefficients: list[float]) -> list[float]:
    """Isolate real roots on [0,1] through derivative roots and bisection."""
    coefficients = coefficients[:]
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    magnitude = max((abs(c) for c in coefficients), default=0)
    if len(coefficients) <= 1 or magnitude == 0:
        return []
    coefficients = [c / magnitude for c in coefficients]
    derivative = [i * coefficients[i] for i in range(1, len(coefficients))]
    boundaries = [0.0] + polynomial_roots(derivative) + [1.0]
    roots = [x for x in boundaries if abs(polynomial_value(coefficients, x)) < 1e-12]
    for a, b in zip(boundaries, boundaries[1:]):
        fa, fb = polynomial_value(coefficients, a), polynomial_value(coefficients, b)
        if abs(fa) < 1e-12 or abs(fb) < 1e-12 or fa * fb >= 0:
            continue
        for _ in range(60):
            mid = (a + b) / 2
            fm = polynomial_value(coefficients, mid)
            if fa * fm <= 0:
                b = mid
            else:
                a, fa = mid, fm
        roots.append((a + b) / 2)
    unique: list[float] = []
    for root in sorted(roots):
        if not unique or abs(root - unique[-1]) > 1e-10:
            unique.append(root)
    return unique


@dataclass
class BeamResult:
    length: float
    ei: float
    supports: list[Support]
    loads: list[Load]
    nodes: list[float]
    displacement: list[float]
    reactions: list[Reaction]
    forces: list[tuple[float, float]]
    couples: list[tuple[float, float]]
    equilibrium_error: tuple[float, float]

    def check_x(self, x: float) -> None:
        finite_number(x, "x")
        if not 0 <= x <= self.length:
            raise ValueError("x must be within the beam.")

    def shear_moment(self, x: float, side: str = "right") -> tuple[float, float]:
        self.check_x(x)
        if side not in ("left", "right"):
            raise ValueError("Side must be left or right.")
        includes = lambda a: a < x or (a == x and side == "right")
        shear = math.fsum(f for a, f in self.forces if includes(a))
        moment = math.fsum(f * (x-a) for a, f in self.forces if includes(a))
        moment -= math.fsum(c for a, c in self.couples if includes(a))
        for load in self.loads:
            if not distributed(load):
                continue
            h = load.b - load.a
            u = min(max(x-load.a, 0.0), h)
            slope = (intensity(load, load.b) - load.value) / h
            shear -= load.value*u + slope*u*u/2
            moment -= load.value*((x-load.a)*u-u*u/2)
            moment -= slope*((x-load.a)*u*u/2-u*u*u/3)
        return shear, moment

    def interval_load(self, index: int) -> tuple[float, float]:
        a, b = self.nodes[index:index+2]
        mid = a + (b-a)/2
        active = [load for load in self.loads if distributed(load) and load.a < mid < load.b]
        return (math.fsum(intensity(load, a) for load in active),
                math.fsum(intensity(load, b) for load in active))

    def deflection(self, x: float) -> tuple[float, float]:
        self.check_x(x)
        index = min(bisect_right(self.nodes, x)-1, len(self.nodes)-2)
        return self.element_field(index, x-self.nodes[index])

    def element_field(self, index: int, t: float) -> tuple[float, float]:
        a, b = self.nodes[index:index+2]
        q0, q1 = self.interval_load(index)
        slope = (q1-q0)/(b-a)
        v0, m0 = self.shear_moment(a)
        y0, theta0 = self.displacement[2*index:2*index+2]
        y = y0 + theta0*t + (m0*t*t/2 + v0*t**3/6 - q0*t**4/24 - slope*t**5/120)/self.ei
        theta = theta0 + (m0*t + v0*t*t/2 - q0*t**3/6 - slope*t**4/24)/self.ei
        return y, theta

    def extrema(self, quantity: str) -> tuple[tuple[float, float], tuple[float, float]]:
        """Return (minimum value,x), (maximum value,x), including jumps."""
        candidates: list[tuple[float, float]] = []
        for x in self.nodes:
            for side in ("left", "right"):
                if quantity == "y":
                    value = self.deflection(x)[0]
                else:
                    value = self.shear_moment(x, side)[0 if quantity == "V" else 1]
                candidates.append((value, x))
        for index, (a, b) in enumerate(zip(self.nodes, self.nodes[1:])):
            h = b-a
            q0, q1 = self.interval_load(index)
            v0, m0 = self.shear_moment(a)
            if quantity == "V":
                coefficients = [q0, q1-q0]
            elif quantity == "M":
                coefficients = [v0, -q0*h, -(q1-q0)*h/2]
            elif quantity == "y":
                theta0 = self.displacement[2*index+1]
                coefficients = [theta0, m0*h/self.ei, v0*h*h/(2*self.ei),
                                -q0*h**3/(6*self.ei), -(q1-q0)*h**3/(24*self.ei)]
            else:
                raise ValueError("Unknown result quantity.")
            for root in polynomial_roots(coefficients):
                x = a + root*h
                value = self.deflection(x)[0] if quantity == "y" else self.shear_moment(x)[0 if quantity == "V" else 1]
                candidates.append((value, x))
        return min(candidates), max(candidates)

    def solution_text(self) -> str:
        lines = [APP_TITLE + " · WORKED SOLUTION", "",
                 f"L = {self.length:g} m     EI = {self.ei:g} kN·m²",
                 "Downward loads positive; upward reactions positive; couples positive CCW.",
                 "Positive M = sagging; y positive upward (negative = downward).", "",
                 "1. MODEL AND LOAD RESULTANTS"]
        for load in self.loads:
            if load.kind == "Moment":
                lines.append(f"Moment at {load.a:g} m: C = {load.value:g} kN·m")
            elif distributed(load):
                force, first = load_resultant(load)
                lines.append(f"{load.kind} [{load.a:g}, {load.b:g}] m: q₀ = {load.value:g}, q₁ = {intensity(load, load.b):g} kN/m")
                lines.append(f"  W = h(q₀+q₁)/2 = {force:.8g} kN; ∫xq dx = {first:.8g} kN·m")
            else:
                lines.append(f"Point at {load.a:g} m: P = {load.value:g} kN")
        lines += ["", "2. DIRECT STIFFNESS", "Node coordinates (m): " + ", ".join(f"{x:g}" for x in self.nodes),
                  "At each node: d = [y, θ]. Pin/roller: y=0. Fixed: y=θ=0.",
                  "For each element of length h:",
                  "k = EI/h³ × [[12,6h,-12,6h], [6h,4h²,-6h,2h²],",
                  "             [-12,-6h,12,-6h], [6h,2h²,-6h,4h²]]",
                  "For linear downward q₀→q₁, equivalent nodal load:",
                  "fₑ = [-h(7q₀+3q₁)/20, -h²(3q₀+2q₁)/60,",
                  "      -h(3q₀+7q₁)/20, +h²(2q₀+3q₁)/60]",
                  "Assemble K and f; solve Kff df = ff; reactions r = Kd − f.",
                  "Diagonal scaling and partial pivoting are used.", "", "3. REACTIONS"]
        for reaction in self.reactions:
            lines.append(f"{reaction.support.kind} at {reaction.support.x:g} m: R = {reaction.force:.8g} kN; C = {reaction.moment:.8g} kN·m")
        downward = math.fsum(load_resultant(load)[0] for load in self.loads)
        first = math.fsum(load_resultant(load)[1] for load in self.loads)
        applied = math.fsum(load.value for load in self.loads if load.kind == "Moment")
        lines += ["", "4. EQUILIBRIUM CHECK",
                  f"ΣR − ΣW = {math.fsum(r.force for r in self.reactions):.8g} − {downward:.8g} = {self.equilibrium_error[0]:.3g} kN",
                  f"Σ(Rx+C) + ΣCload − ∫xq dx − ΣPa = {first-applied+self.equilibrium_error[1]:.8g} + {applied:.8g} − {first:.8g} = {self.equilibrium_error[1]:.3g} kN·m",
                  "", "5. SHEAR AND MOMENT", "H is the step function; ⟨z⟩ = max(z,0). F is signed upward.",
                  "V(x) = ΣF H(x−a) − ∫₀ˣq(t)dt",
                  "M(x) = ΣF⟨x−a⟩ − ΣC H(x−a) − ∫₀ˣq(t)(x−t)dt",
                  "A downward point load jumps V by −P; a CCW couple jumps M by −C.",
                  "Between breakpoints: V′ = −q and M′ = V.",
                  "Extrema include both sides of breakpoints and roots of q=0 / V=0."]
        for quantity, unit in (("V", "kN"), ("M", "kN·m"), ("y", "m")):
            low, high = self.extrema(quantity)
            lines.append(f"{quantity}: min {low[0]:.8g} {unit} at {low[1]:.8g} m; max {high[0]:.8g} {unit} at {high[1]:.8g} m")
        lines += ["", "6. DEFLECTION", "EI y″ = M. Inside each element, t=x−x₀ and s=(q₁−q₀)/h:",
                  "y = y₀ + θ₀t + (M₀t²/2 + V₀t³/6 − q₀t⁴/24 − st⁵/120)/EI",
                  "θ = θ₀ + (M₀t + V₀t²/2 − q₀t³/6 − st⁴/24)/EI",
                  "Deflection extrema solve θ=0 using polynomial root isolation.", "",
                  "Node          y (m)          θ (rad)"]
        for i, x in enumerate(self.nodes):
            lines.append(f"{x:<12.6g} {self.displacement[2*i]:<14.7g} {self.displacement[2*i+1]:.7g}")
        lines += ["", "Assumptions: constant EI, small elastic deflection, ideal supports, bending only.",
                  "Negative reactions may indicate uplift; contact loss is not modelled.",
                  "Educational analysis only; no member strength or building-code checks."]
        return "\n".join(lines)


def solve_beam(length: float, ei: float, supports: list[Support],
               loads: list[Load]) -> BeamResult:
    validate_model(length, ei, supports, loads)
    nodes = sorted({0.0, length, *(s.x for s in supports), *(load.a for load in loads),
                    *(load.b for load in loads if distributed(load))})
    index = {x: i for i, x in enumerate(nodes)}
    size = 2 * len(nodes)
    stiffness = [[0.0] * size for _ in range(size)]
    force = [0.0] * size
    try:
        for i, (a, b) in enumerate(zip(nodes, nodes[1:])):
            h = b-a
            ke = element_stiffness(h, ei)
            mid = a+h/2
            active = [load for load in loads if distributed(load) and load.a < mid < load.b]
            q0 = math.fsum(intensity(load, a) for load in active)
            q1 = math.fsum(intensity(load, b) for load in active)
            fe = consistent_load(h, q0, q1)
            dofs = [2*i, 2*i+1, 2*i+2, 2*i+3]
            for local, global_dof in enumerate(dofs):
                force[global_dof] += fe[local]
                for local_j, global_j in enumerate(dofs):
                    stiffness[global_dof][global_j] += ke[local][local_j]
        for load in loads:
            if load.kind == "Point":
                force[2*index[load.a]] -= load.value
            elif load.kind == "Moment":
                force[2*index[load.a]+1] += load.value
        if any(not math.isfinite(v) for row in stiffness for v in row) or any(not math.isfinite(v) for v in force):
            raise ValueError("The model exceeds floating-point range; check lengths, EI and loads.")
        constrained = {2*index[s.x] for s in supports}
        constrained.update(2*index[s.x]+1 for s in supports if s.kind == "Fixed")
        free = [i for i in range(size) if i not in constrained]
        answer = dense_solve([[stiffness[i][j] for j in free] for i in free], [force[i] for i in free])
        displacement = [0.0] * size
        for i, value in zip(free, answer):
            displacement[i] = value
        residual = [math.fsum(stiffness[i][j]*displacement[j] for j in range(size))-force[i]
                    for i in range(size)]
        reactions = [Reaction(s, residual[2*index[s.x]], residual[2*index[s.x]+1] if s.kind == "Fixed" else 0.0)
                     for s in sorted(supports, key=lambda support: support.x)]
        actions = [(r.support.x, r.force) for r in reactions]
        actions.extend((load.a, -load.value) for load in loads if load.kind == "Point")
        couples = [(r.support.x, r.moment) for r in reactions if r.support.kind == "Fixed"]
        couples.extend((load.a, load.value) for load in loads if load.kind == "Moment")
        # Point actions are already in actions, so only distributed forces remain to subtract.
        sum_force = math.fsum(f for _, f in actions) - math.fsum(load_resultant(load)[0] for load in loads if distributed(load))
        sum_moment = math.fsum(x*f for x, f in actions) + math.fsum(c for _, c in couples)
        sum_moment -= math.fsum(load_resultant(load)[1] for load in loads if distributed(load))
        # Absolute intensity bounds remain meaningful when signed loads cancel.
        load_bounds = [(load, (load.b-load.a)*(abs(load.value)+abs(intensity(load, load.b)))/2)
                       for load in loads if distributed(load)]
        force_bound = math.fsum(abs(f) for _, f in actions) + math.fsum(bound for _, bound in load_bounds)
        moment_bound = math.fsum(abs(x*f) for x, f in actions) + math.fsum(abs(c) for _, c in couples)
        moment_bound += math.fsum(load.b*bound for load, bound in load_bounds)
        if abs(sum_force) > 1e-7*max(force_bound, 1e-300) or abs(sum_moment) > 1e-7*max(moment_bound, 1e-300):
            raise ValueError("The solution failed equilibrium; this model is numerically unreliable.")
        result = BeamResult(length, ei, list(supports), list(loads), nodes, displacement,
                            reactions, actions, couples, (sum_force, sum_moment))
        # Integrated fields must agree with the solved nodal displacement and rotation.
        for i, (a, b) in enumerate(zip(nodes, nodes[1:])):
            h = b-a
            y, theta = result.element_field(i, h)
            v0, m0 = result.shear_moment(a)
            q0, q1 = result.interval_load(i)
            slope = (q1-q0)/h
            y_scale = abs(displacement[2*i]) + abs(displacement[2*i+1]*h)
            y_scale += (abs(m0)*h*h/2 + abs(v0)*h**3/6 + abs(q0)*h**4/24 + abs(slope)*h**5/120)/ei
            theta_scale = abs(displacement[2*i+1])
            theta_scale += (abs(m0)*h + abs(v0)*h*h/2 + abs(q0)*h**3/6 + abs(slope)*h**4/24)/ei
            for value, expected, scale in ((y, displacement[2*i+2], y_scale), (theta, displacement[2*i+3], theta_scale)):
                if not math.isfinite(value) or abs(value-expected) > 1e-7*max(scale, abs(expected), 1e-300):
                    raise ValueError("Displacement recovery failed; the model is numerically unreliable.")
        return result
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("The numerical scale is too extreme; check lengths, EI and loads.") from exc


class BeamApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.supports: list[Support] = []
        self.loads: list[Load] = []
        self.result: BeamResult | None = None
        self.solution_window: tk.Toplevel | None = None
        root.title(APP_TITLE)
        root.geometry("1240x880")
        root.minsize(1020, 760)
        style = ttk.Style(root)
        style.theme_use("clam")
        style.configure("TFrame", background="#f5f6f8")
        style.configure("TLabel", background="#f5f6f8", foreground="#253348")
        style.configure("TLabelframe", background="#f5f6f8")
        style.configure("TLabelframe.Label", background="#f5f6f8", foreground="#253348")
        style.configure("TButton", padding=(9, 5))
        style.configure("Treeview", rowheight=25)
        root.configure(background="#f5f6f8")
        self.length = tk.StringVar(value="6")
        self.ei = tk.StringVar(value="30000")
        self.stiffness_mode = tk.StringVar(value="Constant EI")
        self.modulus = tk.StringVar(value="200")
        self.inertia = tk.StringVar(value="150000000")
        self.support_kind = tk.StringVar(value="Pinned")
        self.support_x = tk.StringVar(value="0")
        self.load_kind = tk.StringVar(value="Point")
        self.load_a = tk.StringVar(value="3")
        self.load_b = tk.StringVar(value="6")
        self.load_value = tk.StringVar(value="30")
        self.load_end = tk.StringVar(value="30")
        self.read_x = tk.StringVar(value="3")
        self.status = tk.StringVar(value="Build a model, then Calculate.")
        self.summary = tk.StringVar(value="No calculation yet.")
        self.readout = tk.StringVar(value="")
        self._build_ui()
        self.length.trace_add("write", self._input_changed)
        self.ei.trace_add("write", self._input_changed)
        self.modulus.trace_add("write", self._input_changed)
        self.inertia.trace_add("write", self._input_changed)
        self.stiffness_mode.trace_add("write", self._stiffness_fields)
        self.load_kind.trace_add("write", self._load_fields)
        self._load_fields()
        self._stiffness_fields()
        self.example()

    def _build_ui(self) -> None:
        main = ttk.Frame(self.root, padding=16)
        main.pack(fill="both", expand=True)
        header = ttk.Frame(main)
        header.pack(fill="x", pady=(0, 12))
        ttk.Button(header, text="About", command=self.show_about).pack(side="right", padx=(12, 0))
        heading = ttk.Frame(header)
        heading.pack(side="left", fill="x", expand=True)
        title_label = ttk.Label(heading, text=APP_TITLE, font=("Segoe UI", 17, "bold"), wraplength=1080)
        title_label.pack(anchor="w", fill="x")
        ttk.Label(heading, text=NAME, font=("Segoe UI", 11)).pack(anchor="w", pady=(3, 0))
        header.bind("<Configure>", lambda event: title_label.configure(wraplength=max(event.width-95, 300)))
        footer = ttk.Frame(main)
        footer.pack(side="bottom", fill="x", pady=(10, 0))
        status_label = ttk.Label(footer, textvariable=self.status, wraplength=880)
        status_label.pack(side="left", fill="x", expand=True)
        footer.bind("<Configure>", lambda event: status_label.configure(wraplength=max(event.width, 200)))
        body = ttk.Frame(main)
        body.pack(fill="both", expand=True)
        setup = ttk.Frame(body, width=360)
        setup.pack(side="left", fill="y", padx=(0, 16))
        setup.pack_propagate(False)
        output = ttk.Frame(body)
        output.pack(side="left", fill="both", expand=True)
        beam = ttk.LabelFrame(setup, text="Beam", padding=10)
        beam.pack(fill="x", pady=(0, 10))
        self._field(beam, "Length (m)", self.length, 0)
        ttk.Label(beam, text="Stiffness").grid(row=1, column=0, sticky="w", pady=3)
        ttk.Combobox(beam, textvariable=self.stiffness_mode, values=("Constant EI", "E and I"), state="readonly", width=16).grid(row=1, column=1, sticky="ew", padx=(8, 0))
        self.ei_entry = self._field(beam, "EI (kN·m²)", self.ei, 2)
        self.modulus_entry = self._field(beam, "E (GPa)", self.modulus, 3)
        self.inertia_entry = self._field(beam, "I (mm⁴)", self.inertia, 4)
        self.stiffness_rows = {"Constant EI": beam.grid_slaves(row=2),
                               "E and I": beam.grid_slaves(row=3) + beam.grid_slaves(row=4)}
        actions = ttk.Frame(setup)
        actions.pack(fill="x", pady=(0, 6))
        for text, command in (("Example", self.example), ("Calculate", self.calculate), ("Clear", self.clear)):
            ttk.Button(actions, text=text, command=command).pack(side="left", padx=(0, 5))
        self.solution_button = ttk.Button(setup, text="Show Solution", command=self.show_solution, state="disabled")
        self.solution_button.pack(fill="x", pady=(0, 12))
        editor = ttk.Notebook(setup)
        editor.pack(fill="x", pady=(0, 10))
        support = ttk.Frame(editor, padding=10)
        editor.add(support, text="Supports")
        ttk.Label(support, text="Type").grid(row=0, column=0, sticky="w", pady=3)
        ttk.Combobox(support, textvariable=self.support_kind, values=("Pinned", "Roller", "Fixed"), state="readonly", width=15).grid(row=0, column=1, sticky="ew", padx=(8, 0))
        self._field(support, "x (m)", self.support_x, 1)
        self._edit_buttons(support, self.add_support, self.save_support, self.remove_support, 2)
        self.support_table = self._table(support, ("Type", "x (m)"), (150, 90), 3, 3)
        self.support_table.bind("<<TreeviewSelect>>", self.select_support)
        load = ttk.Frame(editor, padding=10)
        editor.add(load, text="Loads")
        ttk.Label(load, text="Type").grid(row=0, column=0, sticky="w")
        ttk.Combobox(load, textvariable=self.load_kind, values=("Point", "Moment", "UDL", "Linear"), state="readonly", width=15).grid(row=0, column=1, sticky="ew", padx=(8, 0))
        self._field(load, "Start / x (m)", self.load_a, 1)
        self.load_b_entry = self._field(load, "End (m)", self.load_b, 2)
        self.value_label = ttk.Label(load)
        self.value_label.grid(row=3, column=0, sticky="w", pady=3)
        ttk.Entry(load, textvariable=self.load_value, width=16).grid(row=3, column=1, sticky="ew", padx=(8, 0))
        self.end_label = ttk.Label(load, text="End q (kN/m)")
        self.end_label.grid(row=4, column=0, sticky="w", pady=3)
        self.end_entry = ttk.Entry(load, textvariable=self.load_end, width=16)
        self.end_entry.grid(row=4, column=1, sticky="ew", padx=(8, 0))
        self._edit_buttons(load, self.add_load, self.save_load, self.remove_load, 5)
        self.load_table = self._table(load, ("Type", "Position (m)", "Value"), (60, 105, 125), 6, 4)
        self.load_table.bind("<<TreeviewSelect>>", self.select_load)
        ttk.Label(setup, text="Forces/q: + downward, − upward\nCouples: + CCW, − clockwise\nLinear: triangle or trapezoid; use start/end q.", wraplength=345).pack(anchor="w", pady=(0, 8))
        self.beam_canvas = tk.Canvas(output, height=155, background="white", highlightthickness=1, highlightbackground="#dce1e8")
        self.beam_canvas.pack(fill="x", pady=(0, 10))
        self.beam_canvas.bind("<Configure>", lambda _event: self.draw_beam())
        reaction_frame = ttk.LabelFrame(output, text="Reactions · + upward / + counterclockwise", padding=8)
        reaction_frame.pack(fill="x", pady=(0, 8))
        self.reaction_table = self._table(reaction_frame, ("Support", "x (m)", "R (kN)", "C (kN·m)"), (100, 80, 125, 125), 0, 3)
        summary_label = ttk.Label(output, textvariable=self.summary, wraplength=790)
        summary_label.pack(fill="x", pady=(0, 5))
        output.bind("<Configure>", lambda event: summary_label.configure(wraplength=max(event.width, 200)))
        read = ttk.Frame(output)
        read.pack(side="bottom", fill="x")
        ttk.Label(read, text="x (m)").pack(side="left")
        ttk.Entry(read, textvariable=self.read_x, width=9).pack(side="left", padx=6)
        self.read_button = ttk.Button(read, text="Read", command=self.read_values, state="disabled")
        self.read_button.pack(side="left")
        read_label = ttk.Label(read, textvariable=self.readout, wraplength=590)
        read_label.pack(side="left", padx=10)
        read.bind("<Configure>", lambda event: read_label.configure(wraplength=max(event.width-190, 200)))
        charts = ttk.Frame(output)
        charts.pack(fill="both", expand=True)
        charts.columnconfigure(0, weight=1)
        self.chart_canvases: dict[str, tk.Canvas] = {}
        for row, quantity in enumerate(("V", "M", "y")):
            charts.rowconfigure(row, weight=1)
            canvas = tk.Canvas(charts, height=110, background="white", highlightthickness=1, highlightbackground="#dce1e8")
            canvas.grid(row=row, column=0, sticky="nsew", pady=(0, 7))
            canvas.bind("<Configure>", lambda _event, q=quantity: self.draw_chart(q))
            self.chart_canvases[quantity] = canvas

    @staticmethod
    def _field(parent: ttk.Frame, label: str, variable: tk.StringVar, row: int) -> ttk.Entry:
        parent.columnconfigure(1, weight=1)
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=3)
        entry = ttk.Entry(parent, textvariable=variable, width=16)
        entry.grid(row=row, column=1, sticky="ew", padx=(8, 0))
        return entry

    @staticmethod
    def _edit_buttons(parent: ttk.Frame, add, save, remove, row: int) -> None:
        buttons = ttk.Frame(parent)
        buttons.grid(row=row, column=0, columnspan=2, sticky="ew", pady=7)
        for text, command in (("Add", add), ("Save", save), ("Remove", remove)):
            ttk.Button(buttons, text=text, command=command).pack(side="left", padx=(0, 4))

    @staticmethod
    def _table(parent: ttk.Frame, columns: tuple[str, ...], widths: tuple[int, ...], row: int, height: int) -> ttk.Treeview:
        frame = ttk.Frame(parent)
        frame.grid(row=row, column=0, columnspan=2, sticky="ew")
        parent.columnconfigure(1, weight=1)
        table = ttk.Treeview(frame, columns=columns, show="headings", height=height, selectmode="browse")
        for column, width in zip(columns, widths):
            table.heading(column, text=column)
            table.column(column, width=width, minwidth=40, anchor="center", stretch=True)
        scroll = ttk.Scrollbar(frame, orient="vertical", command=table.yview)
        table.configure(yscrollcommand=scroll.set)
        table.pack(side="left", fill="x", expand=True)
        scroll.pack(side="right", fill="y")
        return table

    def _load_fields(self, *_args) -> None:
        kind = self.load_kind.get()
        is_distributed = kind in ("UDL", "Linear")
        self.load_b_entry.configure(state="normal" if is_distributed else "disabled")
        self.end_entry.configure(state="normal" if kind == "Linear" else "disabled")
        self.value_label.configure(text={"Point": "P (kN)", "Moment": "C (kN·m)", "UDL": "q (kN/m)", "Linear": "Start q (kN/m)"}[kind])

    def _stiffness_fields(self, *_args) -> None:
        mode = self.stiffness_mode.get()
        for name, widgets in self.stiffness_rows.items():
            for widget in widgets:
                if name == mode:
                    widget.grid()
                else:
                    widget.grid_remove()
        self._input_changed()

    def _input_changed(self, *_args) -> None:
        self.invalidate()
        self.draw_beam()

    def invalidate(self) -> None:
        self.result = None
        self.solution_button.configure(state="disabled")
        self.read_button.configure(state="disabled")
        self.reaction_table.delete(*self.reaction_table.get_children())
        self.summary.set("Model changed. Calculate to update results.")
        self.readout.set("")
        self.status.set("Model changed. Calculate to update results.")
        if self.solution_window and self.solution_window.winfo_exists():
            self.solution_window.destroy()
        self.solution_window = None
        for quantity in self.chart_canvases:
            self.draw_chart(quantity)
        self.draw_beam()

    def refresh_model(self) -> None:
        self.support_table.delete(*self.support_table.get_children())
        for i, support in enumerate(self.supports):
            self.support_table.insert("", "end", iid=str(i), values=(support.kind, f"{support.x:g}"))
        self.load_table.delete(*self.load_table.get_children())
        for i, load in enumerate(self.loads):
            position = f"{load.a:g}–{load.b:g}" if distributed(load) else f"{load.a:g}"
            value = f"{load.value:g}→{load.end_value:g}" if load.kind == "Linear" else f"{load.value:g}"
            self.load_table.insert("", "end", iid=str(i), values=(load.kind, position, value))
        self.invalidate()
        self.draw_beam()

    def beam_parameters(self) -> tuple[float, float]:
        return (finite_number(self.length.get(), "Length"),
                stiffness_value(self.stiffness_mode.get(), self.ei.get(), self.modulus.get(), self.inertia.get()))

    def support_input(self) -> Support:
        return Support(self.support_kind.get(), finite_number(self.support_x.get(), "Support x"))

    def load_input(self) -> Load:
        kind = self.load_kind.get()
        value = finite_number(self.load_value.get(), "Load value")
        return Load(kind, finite_number(self.load_a.get(), "Load x/start"), value,
                    finite_number(self.load_b.get(), "Load end") if kind in ("UDL", "Linear") else 0,
                    finite_number(self.load_end.get(), "End intensity") if kind == "Linear" else value)

    def _check_edit(self, supports: list[Support], loads: list[Load]) -> None:
        length, ei = self.beam_parameters()
        # Permit an incomplete support model during editing, but validate every row.
        temporary = supports[:]
        if not any(s.kind == "Fixed" for s in temporary) and len(temporary) < 2:
            used = {s.x for s in temporary}
            for x in (0.0, length):
                if x not in used:
                    temporary.append(Support("Roller", x))
                    used.add(x)
        validate_model(length, ei, temporary, loads)

    def add_support(self) -> None:
        try:
            supports = self.supports + [self.support_input()]
            self._check_edit(supports, self.loads)
            self.supports = supports
            self.refresh_model()
        except ValueError as exc:
            self.error(exc)

    def save_support(self) -> None:
        selection = self.support_table.selection()
        if not selection:
            self.error(ValueError("Select a support row before Save."))
            return
        try:
            supports = self.supports[:]
            supports[int(selection[0])] = self.support_input()
            self._check_edit(supports, self.loads)
            self.supports = supports
            self.refresh_model()
        except ValueError as exc:
            self.error(exc)

    def remove_support(self) -> None:
        selection = self.support_table.selection()
        if selection:
            del self.supports[int(selection[0])]
            self.refresh_model()
        else:
            self.error(ValueError("Select a support row before Remove."))

    def select_support(self, _event) -> None:
        selection = self.support_table.selection()
        if selection:
            support = self.supports[int(selection[0])]
            self.support_kind.set(support.kind)
            self.support_x.set(f"{support.x:g}")

    def add_load(self) -> None:
        try:
            loads = self.loads + [self.load_input()]
            self._check_edit(self.supports, loads)
            self.loads = loads
            self.refresh_model()
        except ValueError as exc:
            self.error(exc)

    def save_load(self) -> None:
        selection = self.load_table.selection()
        if not selection:
            self.error(ValueError("Select a load row before Save."))
            return
        try:
            loads = self.loads[:]
            loads[int(selection[0])] = self.load_input()
            self._check_edit(self.supports, loads)
            self.loads = loads
            self.refresh_model()
        except ValueError as exc:
            self.error(exc)

    def remove_load(self) -> None:
        selection = self.load_table.selection()
        if selection:
            del self.loads[int(selection[0])]
            self.refresh_model()
        else:
            self.error(ValueError("Select a load row before Remove."))

    def select_load(self, _event) -> None:
        selection = self.load_table.selection()
        if selection:
            load = self.loads[int(selection[0])]
            self.load_kind.set(load.kind)
            self.load_a.set(f"{load.a:g}")
            self.load_b.set(f"{load.b:g}")
            self.load_value.set(f"{load.value:g}")
            self.load_end.set(f"{load.end_value:g}")

    def example(self) -> None:
        self.length.set("6")
        self.stiffness_mode.set("Constant EI")
        self.ei.set("30000")
        self.modulus.set("200")
        self.inertia.set("150000000")
        self.supports = [Support("Pinned", 0), Support("Roller", 6)]
        self.loads = [Load("Point", 3, 30)]
        self.read_x.set("3")
        self.refresh_model()
        self.calculate()

    def clear(self) -> None:
        self.supports.clear()
        self.loads.clear()
        self.refresh_model()
        self.status.set("Empty model. Add supports and loads, then Calculate.")

    def calculate(self) -> None:
        self.invalidate()
        try:
            length, ei = self.beam_parameters()
            self.result = solve_beam(length, ei, self.supports, self.loads)
            for reaction in self.result.reactions:
                self.reaction_table.insert("", "end", values=(reaction.support.kind, f"{reaction.support.x:g}", f"{reaction.force:.6g}", f"{reaction.moment:.6g}"))
            moments = self.result.extrema("M")
            deflections = self.result.extrema("y")
            peak_m = max(moments, key=lambda item: abs(item[0]))
            peak_y = max(deflections, key=lambda item: abs(item[0]))
            self.summary.set(f"Peak M: {peak_m[0]:.6g} kN·m at {peak_m[1]:.6g} m   ·   Peak y: {peak_y[0]*1000:.6g} mm at {peak_y[1]:.6g} m")
            self.solution_button.configure(state="normal")
            self.read_button.configure(state="normal")
            for quantity in self.chart_canvases:
                self.draw_chart(quantity)
            self.draw_beam()
            self.read_values()
            self.status.set("Calculated. Ideal constant-EI elastic beam; +M sagging, +y upward. Show Solution explains the formulas.")
        except ValueError as exc:
            self.invalidate()
            self.error(exc)

    def read_values(self) -> None:
        if not self.result:
            return
        try:
            x = finite_number(self.read_x.get(), "x")
            left = self.result.shear_moment(x, "left")
            right = self.result.shear_moment(x)
            y, theta = self.result.deflection(x)
            def pair(a: float, b: float) -> str:
                return f"{a:.6g} / {b:.6g}" if not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9) else f"{b:.6g}"
            self.readout.set(f"V {pair(left[0], right[0])} kN   M {pair(left[1], right[1])} kN·m   y {y*1000:.6g} mm   θ {theta:.6g} rad")
            self.status.set("Read shows left / right values at jumps. Deflection is positive upward.")
        except ValueError as exc:
            self.error(exc)

    def show_solution(self) -> None:
        if not self.result:
            return
        if self.solution_window and self.solution_window.winfo_exists():
            self.solution_window.lift()
            return
        window = tk.Toplevel(self.root)
        self.solution_window = window
        window.title(APP_TITLE + " · Solution")
        window.geometry("850x700")
        text = tk.Text(window, wrap="word", font=("Consolas", 11), padx=18, pady=16, background="white")
        scrollbar = ttk.Scrollbar(window, command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        text.pack(fill="both", expand=True)
        stiffness_note = ""
        if self.stiffness_mode.get() == "E and I":
            stiffness_note = (f"STIFFNESS INPUT\nEI = E × I × 10⁻⁶\n"
                              f"= {self.modulus.get()} GPa × {self.inertia.get()} mm⁴ × 10⁻⁶"
                              f" = {self.result.ei:.8g} kN·m²\n\n")
        text.insert("1.0", stiffness_note + self.result.solution_text())
        text.configure(state="disabled")

    def show_about(self) -> None:
        messagebox.showinfo(APP_TITLE, APP_DESCRIPTION + "\n\n" + NAME, parent=self.root)

    def error(self, exc: ValueError) -> None:
        self.status.set(str(exc))
        messagebox.showerror("Check model", str(exc), parent=self.root)

    def draw_beam(self) -> None:
        canvas = self.beam_canvas
        canvas.delete("all")
        width = max(canvas.winfo_width(), 400)
        try:
            length = finite_number(self.length.get(), "Length")
            if length <= 0:
                return
        except ValueError:
            return
        left, right, y = 48, width-45, 92
        position = lambda x: left + (right-left)*x/length
        canvas.create_text(12, 10, text="MODEL · m", anchor="nw", fill="#253348", font=("Segoe UI", 10, "bold"))
        canvas.create_line(left, y, right, y, width=5, fill="#253348")
        if self.result:
            extreme = max(abs(pair[0]) for pair in self.result.extrema("y"))
            if extreme > 0:
                points: list[float] = []
                for a, b in zip(self.result.nodes, self.result.nodes[1:]):
                    for step in range(25):
                        coordinate = a + (b-a)*step/24
                        displacement = self.result.deflection(coordinate)[0]
                        points.extend((position(coordinate), y-25*displacement/extreme))
                canvas.create_line(*points, fill="#168377", width=2, dash=(5, 3))
            canvas.create_text(right, 12, text="Deformed shape · magnified", anchor="ne", fill="#168377", font=("Segoe UI", 9))
        canvas.create_text(left, 142, text="0", fill="#536078")
        canvas.create_text(right, 142, text=f"{length:g} m", fill="#536078")
        for support in self.supports:
            if not 0 <= support.x <= length:
                continue
            x = position(support.x)
            if support.kind == "Fixed":
                canvas.create_line(x, y-18, x, y+21, fill="#253348", width=4)
                sign = -1 if support.x == 0 else 1
                for offset in range(-18, 20, 8):
                    canvas.create_line(x, y+offset, x+sign*10, y+offset+8, fill="#253348")
            else:
                canvas.create_polygon(x, y+2, x-12, y+20, x+12, y+20, fill="", outline="#253348", width=2)
                if support.kind == "Roller":
                    for dx in (-7, 7):
                        canvas.create_oval(x+dx-3, y+22, x+dx+3, y+28, outline="#253348")
            canvas.create_text(x, y+39, text=f"{support.kind[0]} {support.x:g}", fill="#536078", font=("Segoe UI", 9))
        max_q = max((max(abs(load.value), abs(intensity(load, load.b))) for load in self.loads if distributed(load)), default=1) or 1
        for load in self.loads:
            if not 0 <= load.a <= length:
                continue
            x = position(load.a)
            color = "#276bba"
            if load.kind == "Point":
                start, end = (y-48, y-5) if load.value >= 0 else (y-5, y-48)
                canvas.create_line(x, start, x, end, arrow="last", width=2, fill=color)
                canvas.create_text(x, y-59, text=f"{load.value:g} kN", fill=color, font=("Segoe UI", 9))
            elif load.kind == "Moment":
                canvas.create_arc(x-15, y-34, x+15, y-4, start=30, extent=280 if load.value >= 0 else -280, style="arc", outline="#a25a21", width=2)
                canvas.create_text(x, y-46, text=f"{load.value:g} kN·m {'↺' if load.value >= 0 else '↻'}", fill="#a25a21", font=("Segoe UI", 9))
            elif load.b <= length:
                polyline: list[float] = []
                for step in range(11):
                    coordinate = load.a + (load.b-load.a)*step/10
                    px = position(coordinate)
                    q = intensity(load, coordinate)
                    py = y-6-(abs(q)/max_q)*38
                    polyline.extend((px, py))
                    if q != 0:
                        start, end = (py, y-5) if q > 0 else (y-5, py)
                        canvas.create_line(px, start, px, end, arrow="last", fill=color)
                canvas.create_line(*polyline, fill=color, width=2)
                label = f"{load.value:g} kN/m" if load.kind == "UDL" else f"{load.value:g}→{load.end_value:g} kN/m"
                canvas.create_text((position(load.a)+position(load.b))/2, y-60, text=label, fill=color, font=("Segoe UI", 9))

    def draw_chart(self, quantity: str) -> None:
        canvas = self.chart_canvases[quantity]
        canvas.delete("all")
        labels = {"V": "SHEAR · kN", "M": "MOMENT · kN·m", "y": "DEFLECTION · mm (+ upward)"}
        canvas.create_text(12, 8, text=labels[quantity], anchor="nw", fill="#253348", font=("Segoe UI", 9, "bold"))
        if not self.result:
            canvas.create_text(20, 48, text="Calculate to show this diagram.", anchor="w", fill="#7b8798")
            return
        result = self.result
        width, height = max(canvas.winfo_width(), 400), max(canvas.winfo_height(), 70)
        left, right, top, bottom = 78, width-23, 32, height-23
        factor = 1000 if quantity == "y" else 1
        minimum, maximum = result.extrema(quantity)
        low, high = min(0, minimum[0]*factor), max(0, maximum[0]*factor)
        span = high-low
        if span == 0:
            low, high, span = -1, 1, 2
        low -= span*.12
        high += span*.12
        px = lambda x: left+(right-left)*x/result.length
        py = lambda value: bottom-(bottom-top)*(value-low)/(high-low)
        canvas.create_line(left, py(0), right, py(0), fill="#aeb8c7")
        canvas.create_text(left-8, py(0), text="0", anchor="e", fill="#536078", font=("Segoe UI", 8))
        for value, x in (minimum, maximum):
            canvas.create_text(left-8, py(value*factor), text=f"{value*factor:.4g}", anchor="e", fill="#536078", font=("Segoe UI", 8))
        points: list[float] = []
        pastel = {"V": "#d8e7fa", "M": "#e4dcf4", "y": "#d7eee9"}[quantity]
        for a, b in zip(result.nodes, result.nodes[1:]):
            interval_points: list[float] = []
            for step in range(25):
                x = a+(b-a)*step/24
                if quantity == "y":
                    value = result.deflection(x)[0]
                else:
                    side = "left" if step == 24 else "right"
                    value = result.shear_moment(x, side)[0 if quantity == "V" else 1]
                interval_points.extend((px(x), py(value*factor)))
            canvas.create_polygon(px(a), py(0), *interval_points, px(b), py(0), fill=pastel, outline="")
            points.extend(interval_points)
        canvas.create_line(left, py(0), right, py(0), fill="#aeb8c7")
        if len(points) >= 4:
            color = {"V": "#276bba", "M": "#7954a6", "y": "#168377"}[quantity]
            canvas.create_line(*points, fill=color, width=2)
            critical = max((minimum, maximum), key=lambda pair: abs(pair[0]))
            canvas.create_oval(px(critical[1])-3, py(critical[0]*factor)-3, px(critical[1])+3, py(critical[0]*factor)+3, fill=color, outline="")
        canvas.create_text(left, height-10, text="0", fill="#536078", font=("Segoe UI", 8))
        canvas.create_text(right, height-10, text=f"{result.length:g} m", fill="#536078", font=("Segoe UI", 8))


def main() -> None:
    root = tk.Tk()
    BeamApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
