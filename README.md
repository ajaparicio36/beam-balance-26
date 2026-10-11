# Digital Multi-Load Beam Structural Analysis and Simulation System

A Python-based program that analyzes a structural beam under different loading conditions, such as point loads, uniformly distributed loads, triangular loads, and applied moments. The program allows users to select different types of supports and automatically calculates the support reactions, shear force, bending moment, slope, and deflection. It also provides a digital visualization of the beam and its deformed shape to help users understand how loads affect a structural beam.

Client: **Emman James P. Faderon**. The name appears in the lower-right footer; **About** shows the program description.

Inspired by the public [Optimal Beam calculator](https://optimalbeam.com/beam-calculator.php). This educational reconstruction uses its own interface and solver. It includes advanced beam and load arrangements without accounts or paid services.

## Run

Install Python **3.11 or newer** with Tkinter (included in the usual Windows Python installer). No pip packages are required.

```powershell
python beam_balance.py
```

With this repository's existing virtual environment:

```powershell
.\venv\Scripts\python.exe beam_balance.py
```

The application is entirely in `beam_balance.py`. The other files are documentation and numerical tests. Importing the application file does not open a window.

## Use

1. Set length in **m**. Choose **Constant EI** to enter flexural rigidity in **kN·m²**, or **E and I** to enter E in **GPa** and I in **mm⁴**. Only the active property inputs are used; both modes represent constant stiffness along the beam.
2. In **Supports**, add Pinned, Roller, or Fixed supports at their x positions. Fixed supports must be at an end. Two distinct pins/rollers or one fixed end restrain the beam.
3. In **Loads**, add Point forces, Moment couples, UDLs, or Linear distributed loads. Linear loads cover triangles and trapezoids through their start/end intensities.
4. Select a table row to edit it, then **Save**, or use **Remove**. **Add** creates a new row.
5. **Calculate** updates reactions, shaded diagrams and the deformed beam overlay. The dashed green deformed shape is **magnified** for visibility. **Read** evaluates x, including slope θ in radians; at jumps the shear/moment values are shown as left / right. **Show Solution** opens the derivation and current model's results.

**Example** loads the 6 m simply supported beam with a 30 kN center force and Constant EI=30,000. The matching E/I defaults are E=200 GPa and I=150,000,000 mm⁴. **Clear** removes supports and loads, keeping beam properties. Changes to the model or stiffness inputs clear old results and the deformed overlay; editing draft support/load fields takes effect when Add/Save is clicked.

Signs: forces and intensities are positive **downward**; couples are positive **counterclockwise**; reactions are positive **upward**. Internal moment is positive **sagging**. Deflection is positive upward, so a downward deflection is negative. Reactions use kN and kN·m; deflection diagrams/readouts use mm. The solution also shows y in m and rotation in radians.

For separate E in GPa and I in mm⁴: `EI (kN·m²) = E × I × 10⁻⁶`. For example, E=200 GPa and I=8,000,000 mm⁴ give EI=1,600 kN·m².

## Worked examples

| Model | Input | Expected result |
|---|---|---|
| Simply supported | L=6, EI=30000; Pinned at 0, Roller at 6; Point 30 at 3 | R=15 kN each; Mmax=45 kN·m at 3; ymid=−4.5 mm |
| Fixed-fixed | L=4, EI=1600; Fixed at 0 and 4; UDL 6 over [0,4] | R=12 kN each; CCW couples +8/−8 kN·m; Mmid=4 kN·m; ymid=−2.5 mm |
| Propped cantilever | L=4; Fixed at 0, Roller at 4; UDL 6 over [0,4] | Rleft=15, Rright=9 kN; Cleft=12 kN·m |
| Continuous beam | L=8; Pinned at 0, Rollers at 4 and 8; UDL 6 over [0,8] | R=9,30,9 kN; M at middle support=−12 kN·m |
| Triangular load | L=6; Pinned at 0, Roller at 6; Linear 0→9 over [0,6] | R=9 and 18 kN; Mmax≈20.7846 kN·m at x≈3.4641 m |

Multiple overlapping loads, overhangs, left/right cantilevers, support-coincident loads, signed upward forces, and applied couples are supported. Pins and rollers have the same vertical restraint in this bending-only model.

## How it calculates

The pure-Python solver places nodes at every support/load boundary, assembles Euler–Bernoulli beam stiffness and exact linear-load vectors, solves free nodal DOFs with scaled partial-pivoted elimination, and recovers reactions from `Kd−f`. Shear and moment come from equilibrium; displacement comes from integrating the exact moment field. Reported extrema use polynomial roots and both sides of discontinuities rather than chart samples.

- [Plan](docs/info/PLAN.md)
- [Program flowchart](docs/info/FLOWCHART.md)
- [Step-by-step program algorithm](docs/info/ALGORITHM.md)

## Verify

```powershell
python -m unittest discover -s tests -v
python -m py_compile beam_balance.py
```

The numerical suite checks closed-form determinate, indeterminate, distributed-load, moment-jump, overlapping-load, sign, and validation cases. Compilation/import checks do not exercise desktop interactions.

## Limits

This is one planar beam with positive constant EI, linear elastic small deflection, and ideal zero-settlement supports. It does not calculate axial force, torsion, section capacity, material yielding, buckling, or building-code compliance. Self-weight must be entered as a load. Negative reactions are mathematically allowed and can indicate uplift; unilateral contact is not modelled.

Dense elimination is suitable for classroom models. Extremely close nodes or extreme input scales can be ill-conditioned; such models are rejected when pivot, residual, equilibrium, or recovery checks fail. Root isolation and all calculations use double-precision floating-point arithmetic.
