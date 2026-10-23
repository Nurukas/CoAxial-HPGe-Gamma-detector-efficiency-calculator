# Theory and derivations

Everything the code does, derived from scratch, so each number in the output can be traced.

## 1. What is computed

The full-energy-peak (FEP) efficiency at energy *E* is the probability that a photon of
energy *E*, emitted by the source, deposits **exactly** *E* in the active germanium:

    ε(E) = P(deposit = E | one photon emitted).

**Monte Carlo estimator.** Each analog history *i* is a Bernoulli trial *Xᵢ ∈ {0, 1}*
(full deposit or not). With *N* histories

    p̂ = (1/N) Σ Xᵢ,     E[p̂] = ε,     Var(p̂) = ε(1 − ε)/N.

The relative standard uncertainty is therefore

    u_r = √((1 − ε)/(ε N))   ⇒   N = (1 − ε)/(ε u_r²).

Example: the glass sample at 1332 keV has ε ≈ 0.009; for u_r = 0.5 %,
N ≈ 0.991/(0.009 · 2.5·10⁻⁵) ≈ 4.4·10⁶ histories. The code reports this binomial
uncertainty (`fep_unc`) and OpenMC's batch-to-batch estimate (`fep_unc_batch`); they agree
(see VALIDATION.md), which also confirms that every history is independent.

Set `simulation.rel_unc_target: 0.005` to let OpenMC run batches until the FEP tally reaches
0.5 % (a tally trigger).

## 2. The pulse-height tally and the full-energy window

OpenMC's `pulse-height` score keeps, for every history and every tallied cell, the sum

    h = Σ_collisions in cell (E_before − E_after) − Σ_photons born in cell E_born
        + (energy of photons killed below the cutoff in the cell).

Electrons and positrons are not transported: their kinetic energy stays where they were
created (with the thick-target bremsstrahlung option the bremsstrahlung photons are created
as secondaries and handled by the second term; annihilation photons likewise).
For a photon that is fully absorbed every term cancels except the source energy, so
*h = E* up to floating-point rounding (∼10⁻⁹ eV).

Three consequences built into the code:

1. **One active cell.** The sum *h* is kept *per cell*. If the active volume were split into
   several cells, a photon that deposits 400 keV in one and 262 keV in another would score
   two partial events instead of one 662 keV event. The active germanium is therefore always
   one cell (its region may be arbitrarily complex).

2. **Bin edges.** OpenMC energy bins are (lower, upper]. A bin edge placed exactly at *E*
   would send each full-energy event to one side or the other depending on the last bit of
   rounding. The code always centres a bin (E − w, E + w] on the line.

3. **The window w is part of the definition of "full energy".** A photon that
   Compton-scatters by a small angle θ before reaching the crystal loses

       ΔE = E − E/(1 + κ(1 − cos θ)) ≈ (E²/m_ec²) · θ²/2,     κ = E/m_ec².

   It is counted as "full energy" if ΔE < w, i.e. for scattering angles up to

       θ_max ≈ √(2 w m_ec²) / E.

   With the default w = 1 eV, θ_max ≈ 1° at 60 keV, where the incoherent scattering function
   suppresses scattering anyway: only photons that reach the crystal *without energy loss*
   count — the definition used by Moens et al., ANGLE, EFFTRAN and GESPECOR. With w = 0.5 keV,
   θ_max ≈ 22° at 60 keV, and small-angle scatter in a thick sample starts to be counted as
   peak. A measured peak area sits between the two, depending on resolution and the fitting
   procedure. VALIDATION.md §4 quantifies it: for the glass sample at 46.5 keV a ±0.5 keV
   window adds 4.1 % to the absolute efficiency and 2.9 % to the transfer factor; at 662 keV
   nothing.

Photon transport is analog (survival biasing off): the pulse-height score adds 1 per
matching history and ignores particle weights, so any weight-based variance reduction would
silently bias it.

## 3. Efficiency transfer — why the ratio is robust

Write the Monte Carlo efficiency as a function of the detector-model parameters θ (dead
layers, crystal position, …): ε^MC(θ̂), with true values θ. The transferred efficiency is

    ε_x = ε_ref^meas · ε_x^MC(θ̂) / ε_ref^MC(θ̂).

Expand to first order in δθ = θ̂ − θ, with g = ∂ ln ε / ∂θ:

    ln ε_x^MC(θ̂) = ln ε_x(θ) + g_x δθ,     ln ε_ref^MC(θ̂) = ln ε_ref(θ) + g_ref δθ,

so

    ln ε_x = ln ε_x(θ) + (g_x − g_ref) δθ.

The error is (g_x − g_ref)·δθ instead of g_x·δθ: everything that affects both geometries the
same way cancels. Moens' factorisation ε = (P/T)·Ω̄/4π makes this concrete: an error in the
intrinsic response (P/T) cancels completely; only the difference in how the two geometries
illuminate the crystal (Ω̄) remains. The measured reference also absorbs effects the model
does not contain at all (real dead-layer profile, charge collection, the peak-fitting
convention) as far as they are geometry independent.

Uncertainties of the three factors are independent:

    u²(ε_x)/ε_x² = u²(ε_ref^meas)/ε_ref^meas² + u²(ε_x^MC)/ε_x^MC² + u²(ε_ref^MC)/ε_ref^MC².

The measured reference curve is fitted by weighted least squares,

    ln ε = Σ_{j=0}^{d} a_j [ln(E/E_c)]^j,   weights (ε/u)²,   E_c = geometric mean energy,

with covariance C = (AᵀWA)⁻¹, inflated by χ²_ν when the scatter exceeds the stated
uncertainties; the fit uncertainty at E is √(xᵀCx). `model_ratio_ref = ε_ref^meas/ε_ref^MC`
is reported as a diagnostic: if it departs from 1 in an energy-dependent way, the model
(usually the dead layer) is off — the transfer still compensates to first order.

## 4. Independent cross-check: the effective solid angle (Moens 1981)

For a straight photon path from T in direction u that crosses the active volume in
intervals [t_k, t_k + ℓ_k], the probability that the *first* (non-coherent) interaction is
in the active germanium is

    P(T, u) = Σ_k exp(−τ(t_k)) [1 − exp(−μ_Ge ℓ_k)],     τ(t) = ∫₀ᵗ μ(s) ds.

(The photon survives to t_k with probability e^{−τ(t_k)} and then interacts inside the
segment with probability 1 − e^{−μℓ}. Moens' eq. (11), F_eff = f₁ + f₂f′, is this sum for a
ray crossing the crystal twice through the core.) The effective solid angle is

    Ω̄ = ⟨ ∫_{4π} P(T, u) dΩ ⟩_T.

`hpgeff.esa` integrates this by Monte Carlo over exact ray–cylinder intersections (no
transport), with μ from the same photon library. Directions are sampled uniformly inside a
cone of half-angle θ_c that contains every straight line from the source to the crystal; a
uniform direction in that cone has density 1/[2π(1 − cos θ_c)] per steradian, hence

    Ω̄/4π = [(1 − cos θ_c)/2] · E_cone[P].

Moens' assumption — the "virtual peak-to-total" ε/(Ω̄/4π) does not depend on geometry —
then predicts the transfer ε_x/ε_ref = Ω̄_x/Ω̄_ref. Agreement of that ratio with the full
Monte Carlo ratio checks positions, units, densities, compositions and source sampling by a
completely different calculation (`hpgeff crosscheck`).

## 5. Triangles ("FED"): what faceting does, in numbers

A CAD export replaces a circle of radius R by n flat facets whose vertices lie on the circle.

**Sagitta** (largest inward displacement of the surface):

    s = R(1 − cos(π/n)) ≈ π²R/(2n²).

**Area/volume deficit:** the inscribed n-gon has area (n/2) R² sin(2π/n), so

    A_n/A = (n/2π) sin(2π/n) ≈ 1 − 2π²/(3n²) = 1 − (4/3)(s/R).

Example from this project's cup: R = 32 mm exported with 0.5 mm deflection gave n ≈ 26,
s = 0.23 mm, and a 0.97 % volume deficit (FreeCAD: 159 288 vs 160 850 mm³) — i.e. 3.8 g of
glass missing at fixed density.

**Estimating s from an STL alone.** Two neighbouring facets on a circle of radius R subtend
α = 2π/n, and their normals differ by the same α. A facet's chord is c = 2R sin(α/2), so
R = c/(2 sin(α/2)) and

    s = R(1 − cos(α/2)) = (c/2) tan(α/4).

The removed volume per facet is the circular segment (≈ ⅔ c s, the parabolic segment) times
the facet length. `mesh_report` evaluates both for every pair of neighbouring triangles whose
dihedral angle is below a crease threshold (30°), using the triangle height above the shared
edge as c. On tessellated cylinders this reproduces the exact deficit to < 5 % (unit test).
It assumes vertices lie on the true surface — true for CAD tessellation, not for 3D scans.

**Mass-conserving density.** With `mass_g`, ρ = m/V_mesh: the object has its real mass even
though its surface moved inwards by up to s. What remains is the small redistribution of
material and emitters, an effect first-order in s/R with a small coefficient — measured in
VALIDATION.md.

**Lost particles** (the main neutronics complaint about CAD) come from leaky meshes and from
neighbouring volumes triangulated differently, not from faceting itself. They are prevented
by refusing non-watertight meshes, by measuring overlaps between meshes, and — best — by not
using triangles where exact surfaces are available (built-in shapes, STEP → CSG).

## 6. Optional angular biasing

If every straight line from any emission point to any point of the crystal lies in a
downward cone of half-angle θ_c, emitting only into that cone and multiplying the result by
f = (1 − cos θ_c)/2 is exact for photons that reach the crystal without scattering — the
only photons that can deposit the full energy, apart from elastic (coherent) scatter of
photons emitted outside the cone back towards the crystal. That contribution is tiny
(coherent scattering is forward-peaked) and is checked in VALIDATION.md. The option is off by
default; it pays off for distant sources.

## 7. Physics limits worth knowing

- Electrons are not transported (OpenMC deposits their energy where they are created, plus
  bremsstrahlung). Electron escape through crystal surfaces is neglected — a sub-percent
  effect on FEP at high energies, and it largely cancels in the transfer.
- The dead layer is a sharp boundary; a real transition zone is absorbed by the measured
  reference through the transfer.
- No true-coincidence summing: the reference curve should come from single-line (or
  TCS-corrected) data, as in your plan.
- Resolution is not simulated; it is not needed with the exact full-energy definition.
