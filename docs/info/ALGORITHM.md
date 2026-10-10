# Beam calculation algorithm

## Model and signs

Length and coordinates are in m; force in kN; distributed intensity in kN/m; moment in kN·m; EI in kN·m². Downward entered forces/intensities are positive. Entered couples and reaction couples are positive counterclockwise. Displacement y is positive upward; rotation θ is positive counterclockwise; internal moment M is positive sagging. Thus EI y″ = M. A negative displayed y denotes downward deflection.

Pinned and roller supports constrain y and permit θ. A fixed support constrains both. In this transverse model pinned and roller restraints have the same mathematical effect; horizontal reactions are outside its scope. At least two distinct vertical restraints, or one fixed restraint, are needed. Redundant supports are allowed. Fixed supports are restricted to beam ends.

## Nodes and stiffness

Place nodes at beam ends, support positions, point loads/couples, and distributed-load boundaries. Every element has two DOFs at each end, `[y1, θ1, y2, θ2]`, and length ℓ. Its matrix is

```text
          [ 12      6ℓ    -12      6ℓ ]
k = EI/ℓ³ [ 6ℓ    4ℓ²    -6ℓ    2ℓ² ]
          [-12     -6ℓ     12     -6ℓ ]
          [ 6ℓ    2ℓ²    -6ℓ    4ℓ² ]
```

For signed downward intensity varying linearly from q₀ to q₁ on an element, its exact upward-DOF equivalent load is

```text
f = [ -ℓ(7q₀+3q₁)/20,
      -ℓ²(3q₀+2q₁)/60,
      -ℓ(3q₀+7q₁)/20,
      +ℓ²(2q₀+3q₁)/60 ]
```

This comes from integrating the four cubic Hermite shape functions against `-q(x)`. Sum intensities for overlapping loads. Add a point force as `-P` at the node's y DOF, and a counterclockwise couple as `+C` at its θ DOF. Assemble global K and f, constrain support DOFs to zero, solve `Kff df = ff`, and recover support reactions from `r = Kd - f`.

Diagonal scaling followed by partial-pivoted Gaussian elimination reduces unit/scale sensitivity. Tiny relative pivots, nonfinite calculations, or excessive residuals are reported as unstable/numerically unreliable. No artificial minimum element length or EI is imposed; floating-point resolution and conditioning still limit extreme combinations.

## Equilibrium diagrams

Let H be a step function and `⟨z⟩ = max(z,0)`. Every vertical action F is signed upward (reactions positive, downward point loads negative), at position a. Every external couple C is signed counterclockwise (including fixed reaction couples).

```text
V(x) = Σ F H(x-a) - ∫₀ˣ q(t) dt
M(x) = Σ F ⟨x-a⟩ - Σ C H(x-a) - ∫₀ˣ q(t)(x-t) dt
```

Consequently a downward point force produces a shear jump `-P`, while a positive counterclockwise couple produces a moment jump `-C`. Between breakpoints, `V′=-q` and `M′=V`. Left/right values are shown at discontinuities.

For a linear load beginning at a, ending at b, `h=b-a`, and `s=(q₁-q₀)/h`, set `u=min(max(x-a,0),h)`. Its shear contribution is `-q₀u - s u²/2`. Its moment contribution is

```text
-q₀((x-a)u-u²/2) - s((x-a)u²/2-u³/3).
```

On each interval, shear is quadratic and moment is cubic. Find shear extrema from `q(x)=0`, and moment extrema from `V(x)=0`, plus both sides of every boundary. Diagram drawing samples these exact functions; reported shear/moment extrema are not sampled estimates.

## Displacement

Within an element use the exact left-end moment and shear, q₀, and load slope s:

```text
EI y(x₀+t) = EI(y₀+θ₀t) + M₀t²/2 + V₀t³/6 - q₀t⁴/24 - s t⁵/120.
```

This integrates `EI y″=M` and preserves the stiffness solution's endpoint y and θ. Unlike cubic interpolation alone, it captures the exact quartic/quintic particular solution under distributed loading. Deflection extrema are found by solving θ=0 on each interval using derivative-root isolation and bisection.

## Checks and limitations

Validate all numbers as finite, `L>0`, `EI>0`, coordinates within `[0,L]`, and nonempty load intervals. Reject duplicate supports at the same position; model upgrades should edit that support. Verify global vertical and moment equilibrium and the free-DOF solver residual. Analytical tests check support conditions, diagram jumps, and closed-form beam cases.

The method assumes constant EI, small deformation, linear elasticity, ideal zero-settlement supports, and bending only. Negative support reactions can indicate uplift; the mathematical support is bilateral and does not simulate contact loss. It is an educational analysis tool, not a structural safety or code compliance check.
