# Beam Balance reconstruction plan

## Goal

Build a beginner-friendly desktop beam calculator inspired by the public [Optimal Beam calculator](https://optimalbeam.com/beam-calculator.php). Keep the application in one Python file and use the installed standard-library Tkinter toolkit. The goal is an educational reconstruction of beam modelling, reactions, and diagrams, not a pixel-for-pixel copy or the reference's commercial services.

## Planned deliverables

1. Write this plan, the flowchart, and the calculation algorithm before implementation.
2. Implement `beam_balance.py`: a Tkinter/ttk form, editable support/load tables, beam sketch, shear and moment diagrams, reaction table, and a readable solution window.
3. Accept pinned, roller, and fixed supports; point forces, applied couples, and uniform/triangular/trapezoidal distributed loads. Supports may form cantilevers, overhangs, propped cantilevers, fixed-fixed beams, or continuous beams.
4. Solve the constant-EI Euler–Bernoulli model using direct stiffness, exact equivalent loads, a scaled dense linear solver, and equilibrium-derived diagrams. Reject unstable or numerically unreliable models.
5. Verify analytical cases with standard-library unit tests, run a compile check, and document installation, examples, signs, limitations, and short button meanings.

## Interface

- Beam length (m) and flexural rigidity EI (kN·m²).
- Support type and position, with Add/Edit/Remove actions.
- Load type; position or interval; signed force, couple, or start/end intensity.
- Short main actions: Example, Calculate, Clear, Show Solution.
- Results: reactions, maximum/minimum moment and shear, diagrams, and values at a requested x. Deflection is included because the stiffness solution already determines it.
- Model changes invalidate displayed calculations so stale results cannot be mistaken for current results.

## Verification and delivery

Numerical tests will cover both cantilevers, simple supports, fixed-fixed and propped beams, two-span continuous beams, triangular/overlapping loads, applied moments, diagram jumps/extrema, and validation. Tests and compilation are non-interactive. Native GUI interaction requires separate confirmation that the shared desktop is available. Configure the requested GitHub origin; do not commit or push.

## Boundaries

One planar prismatic beam, positive constant EI, linear elastic small deflection, ideal supports, no support settlement, no axial/torsional response, no self-weight unless entered as a load, and no section capacity/design checks. Dense elimination is appropriate for a classroom model; ill-conditioned configurations are rejected rather than silently reported.
