# Alternative Track Parametrization for ACTS

## Problem

ACTS hard-codes a single bound track parametrization: `(loc0, loc1, φ, θ, q/p, t)`. The direction encoding `d = (sin θ cos φ, sin θ sin φ, cos θ)` has a coordinate pole at `θ = 0, π` — at those points φ is undefined and the free↔bound Jacobian has a `1/sin θ` term that diverges (`Core/include/Acts/Utilities/JacobianHelpers.hpp`, `freeToSphericalDirectionJacobian`).

This is fine for collider geometries where tracks never sit on the z-axis, but it breaks down for telescope detectors with sensors stacked along z, where tracks point along the singular direction by construction. The current workaround — rotating the geometry — is awkward and doesn't scale.

The test suite acknowledges the issue:

- `Tests/IntegrationTests/PropagationEigenConstant.cpp:115` — *"Covariance transport does not work for theta close to poles"*.
- `Tests/IntegrationTests/PropagationTests.hpp:306` — manually shifts disc surfaces off-axis with a `1_cm` offset to dodge `ρ=0, φ-undefined`.

Goal: support a non-singular parametrization for telescope-like geometries, without templating the entire framework on a parametrization trait, and while allowing mixed geometries (barrel cylinders + forward telescope arms) in a single binary.

## Recommended replacement parametrization

For a planar surface with local frame `(e_u, e_v, n)`, encode the direction as **plane-local slopes**:

```
(loc0, loc1, t_x = (d · e_u) / (d · n), t_y = (d · e_v) / (d · n), q/p, t)
```

Same vector dimension as the spherical set — only the meaning of indices 2 and 3 changes. The only singularity is `d · n = 0` (track grazing the sensor plane), which is a real geometric condition rather than a coordinate artifact. This is the same chart used by EUDAQ, Proteus, EUTelescope, and General Broken Lines for beam-telescope fits.

For a disc the position parametrization stays polar `(r, φ_loc)` — slopes for direction compose cleanly with polar position because the disc's tangent space is still a flat 2D plane with a Cartesian frame from `Surface::transform()`. Position↔Cartesian Jacobian sits in the upper-left 2×2 block; slope↔direction Jacobian sits in the 2×2 block at indices (2,3). They don't interact.

## Design

### Core idea

The choice of chart is a property of the **surface**, dispatched at runtime through a **chart policy object** owned by the surface. The free 8-vector inside the propagator stays parametrization-agnostic, so all framework-level machinery (steppers, KF/CKF/GX2F/GSF, MultiTrajectory, vertexing) sees chart-neutral state; the chart is consulted only at the bound↔free boundary and by a small set of chart-aware helpers.

### Vector layout

`BoundVector` stays 6-D. Indices 2 and 3 keep their names (`eBoundPhi`, `eBoundTheta`) for ABI continuity; semantics are reinterpreted to `(t_x, t_y)` when a surface's chart says so. Add neutral aliases `eBoundDir0`, `eBoundDir1` for slope-aware code paths.

### Chart policy object

A single small interface owns everything chart-specific:

```cpp
class BoundParametrizationChart {
public:
  virtual ~BoundParametrizationChart() = default;

  // direction encoding (indices 2-3 of the bound vector)
  virtual Vector3 toDirection(const BoundVector& p, const Transform3& localToGlobal) const = 0;
  virtual void   fromDirection(const Vector3& d, const Transform3& localToGlobal,
                               BoundVector& p) const = 0;

  // direction block of the bound↔free Jacobian
  virtual Eigen::Matrix<double,3,2> directionJacobian(const Vector3& d,
                                                     const Transform3& l2g) const = 0;
  virtual Eigen::Matrix<double,2,3> inverseDirectionJacobian(const Vector3& d,
                                                            const Transform3& l2g) const = 0;

  // periodicity layer (replaces the loose helpers in TrackParameterHelpers)
  virtual void        normalize(BoundVector& p) const = 0;
  virtual BoundVector subtract(const BoundVector& a, const BoundVector& b) const = 0;

  // per-index predicate, needed by GSF merging and diagnostic writers
  virtual bool isPeriodic(BoundIndices i) const = 0;

  // human-readable names for print/IO
  virtual std::array<std::string_view, 6> parameterNames() const = 0;
};

class SphericalChart : public BoundParametrizationChart { /* current behavior */ };
class PlaneSlopeChart : public BoundParametrizationChart { /* slope math, identity normalize */ };
```

Both implementations are stateless singletons (`SphericalChart::instance()`, `PlaneSlopeChart::instance()`).

### Where the chart lives

Exactly one new virtual on `Surface`:

```cpp
class Surface {
  // existing boundToFreeJacobian / freeToBoundJacobian stay, but their
  // internals delegate the direction block to chart()
  virtual const BoundParametrizationChart& chart() const {
    return SphericalChart::instance();          // default
  }
};

class PlaneSurface : public Surface {
  const BoundParametrizationChart* m_chart;     // set in ctor, immutable
  const BoundParametrizationChart& chart() const override { return *m_chart; }
};

class DiscSurface : public Surface {
  const BoundParametrizationChart* m_chart;
  const BoundParametrizationChart& chart() const override { return *m_chart; }
};

class CurvilinearSurface : public Surface {
  const BoundParametrizationChart* m_chart;     // supplied by PropagatorOptions
  const BoundParametrizationChart& chart() const override { return *m_chart; }
};
```

Cylinder, line, straw, perigee, cone: unchanged — no override, inherit the spherical default. They don't have to "advertise" anything; they just are what they are.

### Curvilinear in slope mode

`CurvilinearSurface` is built at a point along the track with normal `n = T` and local axes `(U, V)` from `createCurvilinearUnitU`. In slope coordinates `(t_x = d·U / d·n, t_y = d·V / d·n)`:

- At construction, slopes are exactly `(0, 0)` because `d` and `n` are parallel by definition.
- The Jacobian is well-conditioned everywhere on the surface — `d·n = 1` at construction and stays near 1 for any small perturbation.

The current pole in the curvilinear representation is in the *direction encoding* (global φ, θ), not in the plane geometry — `createCurvilinearUnitU` already handles the z-axis fallback for the local frame. Slopes fix the pole for the default track representation, which is the deeper source of `"covariance transport does not work for theta close to poles"`.

Because curvilinear surfaces are built implicitly by the propagator/stepper (not by user geometry code), the chart choice is made at the **propagator level**, not at surface construction. Add a `curvilinearChart` option to `PropagatorOptions`, default `eSpherical` for backward compatibility, threaded into every site that builds a curvilinear surface:

- `Core/src/Propagator/detail/JacobianEngine.cpp` — `boundToCurvilinearTransportJacobian`
- `Core/include/Acts/Propagator/RiddersStepper.hpp:581,628` — inline `CurvilinearSurface(position, direction)` construction sites
- **`stepper.transportCovarianceToCurvilinear(...)` on every stepper** — this is the one the plan cannot hand-wave. It is called from `MaterialInteractor.hpp:100` (and elsewhere) and it materializes the curvilinear covariance whose direction block the scattering increment lands on. For the chart choice to reach it, the **stepper state must carry the curvilinear chart** (threaded in from `PropagatorOptions` at stepper initialization), because the stepper has no `PropagatorOptions` at the point of that call. This touches `SympyStepper`, `StraightLineStepper`, `EigenStepper`, `AtlasStepper` state initialization, not just two jacobian engines.
- Seeders that emit curvilinear parameters

Rollout: opt-in via `PropagatorOptions` while the spherical default stays in place, exercise the slope variant across the test matrix and downstream pipelines, then flip the default in a later release.

## Mixed parametrization

A single track can thread through spherical-mode cylinders and slope-mode planes/discs in one propagation and one fit, because each surface-crossing is already a round trip `bound_A → free → bound_B` with each surface providing its own chart's conversion. The Jacobian composition

```
J_AB = (surface B chart inverseDirectionJacobian) · (free transport) · (surface A chart directionJacobian)
```

absorbs the chart difference automatically. No "conversion shim" between charts is needed; the free 8-vector is the common interface.

Specifically:
- **Stepper**: untouched. Carries `FreeVector` state; only sees bound at boundaries.
- **CKF / KF / GX2F / GSF predict step**: bound state at each surface comes out in that surface's chart automatically.
- **Measurement update**: `H` operates on bound indices and is chart-agnostic (see next section).
- **Material effects**: applied in free space (see below).
- **Seeding**: produces curvilinear (spherical) parameters. First real surface intersection converts via that surface's own chart.
- **Vertexing**: inputs are perigee or curvilinear, both spherical. Telescope half of a track is unaffected.

## Measurement projector

The projector is row-selection on the bound vector — pure linear algebra, no semantic inspection. The KF update needs only that the measurement vector, its covariance, and the projection share a chart with the surface's predicted state.

- **Position-only measurements** (loc0/loc1) — the entire silicon telescope use case — are identical in both modes. Slope-vs-spherical doesn't enter; position↔direction covariance correlations carry information regardless of which chart the direction lives in.
- **Direction measurements** (rare: segment-based detectors, cosmics, pseudo-measurements) must be expressed as `(t_x, t_y)` on slope-mode surfaces and `(φ, θ)` on spherical ones. Whatever produces the measurement has the surface pointer and can dispatch.
- **Mixed projections** within one measurement work without modification.

No changes required to the projector encoding, the subspace-indices helper, the calibrator interface, or any downstream KF/CKF machinery.

## Angle periodicity handling

The fitters call a periodicity layer at specific points that would silently mangle slope values if left alone. All of these funnel through `chart.normalize`, `chart.subtract`, and `chart.isPeriodic`:

| Site | File:line | Fix |
|---|---|---|
| KF gain-matrix update, post-update normalize | `Core/include/Acts/TrackFitting/detail/GainMatrixUpdaterImpl.hpp:67` | Call `surface.chart().normalize(filtered)` |
| Gain-matrix smoother, innovation | `Core/src/TrackFitting/GainMatrixSmoother.cpp:66` | Call `surface.chart().subtract(prevSmoothed, prevPredicted)` |
| Gain-matrix smoother, normalize | `Core/src/TrackFitting/GainMatrixSmoother.cpp:68` | Call `surface.chart().normalize(smoothed)` |
| MBF smoother, normalize | `Core/src/TrackFitting/MbfSmoother.cpp:26` | Call `surface.chart().normalize(smoothed)` |
| GSF component merging | `Core/include/Acts/TrackFitting/detail/GsfComponentMerging.hpp:130-207` | Per-index dispatch on `chart.isPeriodic(i)`: polar-mean for periodic, arithmetic-mean for non-periodic |
| RiddersStepper numerical Jacobian | `Core/include/Acts/Propagator/RiddersStepper.hpp:533-535` | Call `chart.subtract` instead of `difference_periodic` |
| Parameter validation | `Core/src/EventData/TrackParameterHelpers.cpp:28-30` | Range checks only when the chart declares the index periodic/bounded |

The free-function helpers `normalizeBoundParameters(BoundVector)` and `subtractBoundParameters(BoundVector, BoundVector)` in `Core/include/Acts/EventData/TrackParameterHelpers.hpp` get surface-aware overloads for callers that have a surface. The surface-less overloads stay as "spherical assumed" for callers that genuinely don't have a surface (mostly tests and pre-fit utilities).

The `chart.isPeriodic(i)` predicate is the right granularity for GSF and validation — it lets per-index code paths stay generic instead of branching on chart type. For `SphericalChart` it returns `true` only for `eBoundPhi`. For `PlaneSlopeChart` it's always `false`. Adding a hypothetical future chart that's periodic somewhere else is one line.

## Material effects — free-space rewrite

The multiple-scattering covariance is currently applied with a chart-specific formula (`Core/src/Propagator/detail/PointwiseMaterialInteraction.cpp:69-75`):

```
varianceTheta = θ_0²
variancePhi   = θ_0² / sin²(θ)
```

and written straight into `cov(eBoundPhi/eBoundTheta)` (`PointwiseMaterialInteraction.hpp:213-218`). Note that in the propagator path this `cov` is the **curvilinear** covariance, not a surface-bound one — `MaterialInteractor.hpp:100` transports to curvilinear first. Either way the formula bakes the spherical chart into the physics model. The fix is to apply scattering in **free space** — where it physically belongs — and let the chart's Jacobian project it into whichever representation the covariance currently lives in.

### Design

Scattering is a rank-2 Gaussian on the plane perpendicular to the track direction. Compute the 3×3 free direction covariance increment `Σ_dir` once (chart-neutral), then transport into bound via a Jacobian sandwich:

```
ΔΣ_bound[2:4, 2:4] = J_dir · Σ_dir · J_dirᵀ
```

where `J_dir = chart.inverseDirectionJacobian(d, l2g)` is the same 2×3 (`d(dir)/d(bound)`) block the chart already provides. Dimensions: `(2×3)·(3×3)·(3×2) = 2×2`, dropped into the direction block of the bound covariance. No inverse, no per-chart math.

> **Note the transpose order.** The projection is `J · Σ · Jᵀ` with `J` the 2×3 free→bound direction Jacobian — *not* `Jᵀ · Σ · J`, which is non-conforming for a 2×3 `J`. Get this backwards and the code either won't compile or (if `J` is swapped for the 3×2 forward Jacobian) produces a 3×3 result that silently corrupts the covariance.

**Which chart applies depends on the call site — this is the subtle part.** In the propagator path, `MaterialInteractor` calls `stepper.transportCovarianceToCurvilinear(state.stepping)` *before* applying the increment (`MaterialInteractor.hpp:100`), so the covariance is in the **curvilinear** representation when scattering lands. The relevant chart there is the *curvilinear* chart (`PropagatorOptions::curvilinearChart`), **not** the surface's chart. In the KF/GSF/GX2F sites the covariance is in the surface-bound representation, so the *surface's* chart applies. The free-space rewrite is correct for both, but the sandwich helper must be handed the right chart per site — do not assume `surface.chart()` everywhere.

`computeMaterialEffects` changes from returning `{var_φ, var_θ}` to returning `Σ_dir`. A new helper (either on the chart or as a free function taking a surface reference) does the sandwich into the bound covariance.

### Blast radius

Five call sites:

- `Core/include/Acts/Propagator/MaterialInteractor.hpp:99-102` — the actor everything else routes through.
- `Core/include/Acts/TrackFitting/KalmanFitter.hpp:463, 557, 638` — three sites in KF.
- `Core/include/Acts/TrackFitting/detail/GsfUtils.hpp:389-404` — GSF component-level scattering.

Each replaces "add `var_φ`/`var_θ` to the bound cov" with one helper call. No data-flow rewiring beyond the signature and the call.

**GX2F simplifies.** `Core/include/Acts/TrackFitting/GlobalChiSquareFitter.hpp:469,481-485` currently multiplies the φ contribution by `sin²θ` precisely because the bound formulation is non-orthonormal in the angle basis. In free space the correction disappears — the code deletes rather than dispatches.

**Plugins unaffected.** EDM4hep, DD4hep, Geant4, FATRAS all provide material *slabs* through `ISurfaceMaterial`; none call `computeMaterialEffects` directly.

### Why not the chart-aware bound formula instead

Alternative: keep scattering in bound space, have each chart supply its own scattering covariance formula. Smaller signature diff, but:

- Scattering is *physically* free-space; the bound formulation is a coordinate accident that only worked when there was one chart.
- Each new chart requires a fresh formula derivation. Recurring correctness exposure.
- GX2F's `sin²θ` math has to stay and become chart-aware, whereas free-space just deletes it.

The free-space rewrite pays the cost once and covers every chart, including future ones. Recommended.

## Accessor semantics — the hardest problem, not a footnote

This is the deepest structural issue in the whole change and it needs to be resolved *before* slope mode is enabled anywhere, not deferred to a late "policy decision."

`GenericBoundTrackParameters` today assumes **the world direction is recoverable from the 6 params alone**. The accessors take no `GeometryContext`:

```cpp
double  phi()       const { return m_params[eBoundPhi]; }          // BoundTrackParameters.hpp:242
double  theta()     const { return m_params[eBoundTheta]; }        //                        :245
Vector3 direction() const { /* from m_params[eBoundPhi/Theta] */ } //                        :253
```

and `position(gctx)`, `momentum()`, `fourPosition(gctx)`, `referenceFrame(gctx)` are all built on `direction()`.

In **spherical** mode this is sound: `(φ, θ)` are global angles, intrinsic to the params. In **slope** mode the direction is defined relative to the surface local frame `(e_u, e_v, n)`, which is **alignment/context-dependent**. So the world direction is a function of the surface `Transform3`, which lives behind a `GeometryContext` the accessor does not have. There is no way to make `direction()` chart-correct without either a context or a cached direction.

Two viable resolutions — pick one project-wide, and decide in PR1, not PR8:

- **(A) Cache the resolved world direction in the parameters object.** At construction (which already takes a `GeometryContext`, `BoundTrackParameters.hpp:59`) resolve slopes→direction once and store the `Vector3`. `direction()`/`phi()`/`theta()` return from the cache; they stay context-free and every downstream caller is unaffected. Cost: 3 extra doubles per params object and an invariant that the cache tracks `m_params`. **Recommended** — it preserves the accessor signatures, which is what makes the 242-site audit tractable.
- **(B) Add a `GeometryContext` parameter to `direction()` and everything downstream of it.** Semantically clean but ripples through a large fraction of the EDM and every caller of `direction()`/`momentum()`/`position()`. Rejected unless (A) proves unworkable.

Under either option `phi()`/`theta()` are defined as **the spherical angles of the world direction** (not raw `params[eBoundPhi]`). Callers that genuinely want the raw slope component read `params()[eBoundDir0]` explicitly.

This decision gates the PR ordering: the accessor/caching work (currently PR8) must move ahead of any slope-enabled surface (PR3), because the moment a slope-mode `PlaneSurface` exists, every context-free `direction()` call on its parameters is wrong.

## IO, print, and EDM4hep boundaries

Several output paths hardcode "phi"/"theta" labels or wrap residuals with `difference_periodic`. These need chart-aware treatment:

- `Core/src/EventData/PrintParameters.cpp:29-30` — `makeBoundNames()` returns hardcoded strings. Route through `surface.chart().parameterNames()`.
- `Examples/Io/Root/src/RootTrackStatesWriter.cpp:580,582` — writes `parameters[eBoundPhi]` to a branch named `"phi"`. Either use chart-aware branch names or convert slope-mode direction into spherical before writing (choose one project-wide).
- `Examples/Io/Root/src/RootTrackStatesWriter.cpp:636-637` — residuals use `difference_periodic` on `eBoundPhi`. Route through `chart.subtract`.
- `Plugins/EDM4hep/src/EDM4hepUtil.cpp:164-165,189-190` — explicit `tan(π/2 - θ)` / `sin(θ)` conversions in `jacobianToEdm4hep`. EDM4hep's schema only supports a helix/perigee parametrization; slope-mode states cannot round-trip. **Policy: convert to spherical at the EDM4hep boundary** (lossy but well-defined, since the world direction fully determines the spherical angles). Document explicitly.

## Touch list

| Area | Action | Risk |
|------|--------|------|
| `Definitions/TrackParametrization.hpp` | Add `eBoundDir0`/`eBoundDir1` aliases | Low |
| `Utilities/JacobianHelpers.hpp` | Add slope↔direction Jacobians alongside spherical ones | Low |
| `Utilities/UnitVectors.hpp` | Add `makeDirectionFromSlopes(...)` | Low |
| **New**: `EventData/BoundParametrizationChart.{hpp,cpp}` + `SphericalChart`, `PlaneSlopeChart` | Interface + both implementations, singleton access | Medium |
| `Surfaces/Surface.hpp` | Add one virtual `chart()`, default returns `SphericalChart::instance()`; existing Jacobians route the direction block through `chart()` | Medium |
| `Surfaces/PlaneSurface.{hpp,cpp}` | Store chart in ctor; override `chart()` | Low |
| `Surfaces/DiscSurface.{hpp,cpp}` | Same pattern; keep polar position | Low |
| `Surfaces/CurvilinearSurface.{hpp,cpp}` | Same pattern; chart supplied by `PropagatorOptions` | Medium |
| `Propagator/PropagatorOptions` | Add `curvilinearChart` field, default `eSpherical` | Low |
| `Propagator/detail/JacobianEngine.cpp` | Thread `curvilinearChart` into `boundToCurvilinearTransportJacobian` | Medium |
| `Propagator/RiddersStepper.hpp:181,533-535,581,628` | Thread `curvilinearChart` into inline `CurvilinearSurface(...)` sites; route φ diff through `chart.subtract` | Medium |
| `EventData/TransformationHelpers.cpp` | Route direction part through `surface.chart()` | Medium |
| `EventData/TrackParameterHelpers.{hpp,cpp}` | Surface-aware overloads of `normalize`/`subtract`/`add`; validation only when chart declares periodic/bounded | Low-medium |
| `TrackFitting/detail/GainMatrixUpdaterImpl.hpp:67` | Call `surface.chart().normalize(...)` | Low |
| `TrackFitting/GainMatrixSmoother.cpp:66,68` | Call `chart.subtract` + `chart.normalize` | Low |
| `TrackFitting/MbfSmoother.cpp:26` | Call `chart.normalize` | Low |
| `TrackFitting/detail/GsfComponentMerging.hpp:130-207` | Per-index dispatch on `chart.isPeriodic(i)` | Medium-high |
| `EventData/BoundTrackParameters` accessors | Cache resolved world direction at construction (resolution A); `phi()`/`theta()`/`direction()` read the cache; single mutation funnel re-resolves | **High** — invariant-critical, lands PR2.5 |
| **Material rewrite**: `Propagator/detail/PointwiseMaterialInteraction.{hpp,cpp}` | Return free-space `Σ_dir` (3×3 rank-2) instead of `{var_φ, var_θ}`; move the write from `cov(eBoundPhi/eBoundTheta)` (`.hpp:213-218`) to the `J·Σ·Jᵀ` sandwich | Medium |
| `Propagator/MaterialInteractor.hpp:99-102` | Sandwich `Σ_dir` via **the curvilinear chart** (cov is transported to curvilinear at `:100`), *not* the surface chart | Medium |
| Stepper state (`SympyStepper`, `StraightLineStepper`, `EigenStepper`, `AtlasStepper`) | Carry `curvilinearChart` in stepper state so `transportCovarianceToCurvilinear` and curvilinear hand-offs use the right chart | Medium |
| `TrackFitting/KalmanFitter.hpp:463, 557, 638` | Same helper call at each site | Low |
| `TrackFitting/detail/GsfUtils.hpp:389-404` | Same helper call per GSF component | Low |
| `TrackFitting/GlobalChiSquareFitter.hpp:469,481-485` | Delete `sin²θ` scaling, use the free-space helper | Medium |
| `EventData/PrintParameters.cpp:29-30` | Route names through `chart.parameterNames()` | Low |
| `Examples/Io/Root/src/RootTrackStatesWriter.cpp:580,582,636-637` | Chart-aware branch names; residuals via `chart.subtract` | Low |
| `Plugins/EDM4hep/src/EDM4hepUtil.cpp:164-165,189-190` | Convert to spherical at the EDM4hep boundary; document as lossy | Low |
| `Examples/Detectors/TelescopeDetector` | Constructor flag to opt sensors into slope mode | Low |
| Python bindings | Expose chart enum + ctor arg on `PlaneSurface`/`DiscSurface`; expose `PropagatorOptions::curvilinearChart` | Low |
| Tests | Parametrize existing angle-bound tests on chart; new slope Jacobian unit tests; forward-track propagation without `1_cm` workaround; scattering free→bound helper test | Medium |
| Steppers (Eigen, Sympy, Atlas, StraightLine, GX2F actor) | **No structural changes** — free-vector state carries through | — |
| MultiTrajectory, TrackContainer | **No changes** | — |
| Cylinder, line, straw, cone surfaces | **No changes** | — |
| Vertexing, seeding | **Audit only** — consumers of curvilinear states routed through accessors that recompute from direction | Medium (in curvilinear phase) |

## Phasing

### PR1 — chart interface and slope helpers

Introduce `BoundParametrizationChart`, `SphericalChart`, `PlaneSlopeChart`. Add slope↔direction Jacobians in `JacobianHelpers.hpp` and `makeDirectionFromSlopes` in `UnitVectors.hpp`. Unit tests: chart methods against numerical derivatives, chart round-tripping. No behavior change anywhere else — nothing consumes the new charts yet.

### PR2 — refactor existing math to route through `SphericalChart`

Add `chart()` virtual on `Surface` (default: spherical singleton). Route `TransformationHelpers`, the four surface Jacobians, the periodic helpers, KF/CKF/GX2F/GSF/MBF/Ridders/validation call sites, and Print/IO through `chart()`. Every code path still resolves to `SphericalChart::instance()`; existing tests stay green bit-for-bit. Reviewer's job is to verify the refactor is faithful.

### PR2.5 — direction resolution / accessor hardening (**must precede PR3**)

Implement resolution (A) from *Accessor semantics*: cache the resolved world direction in `GenericBoundTrackParameters`, resolved at construction via `surface.chart().toDirection(...)`, and redefine `phi()`/`theta()`/`direction()` to read the cache. Still spherical-only at this point, so the cache is exactly today's `makeDirectionFromPhiTheta` — no behavior change, fully regression-testable. Add `eBoundDir0`/`eBoundDir1` aliases and do the grep-and-review pass over raw `params[eBoundPhi]`/`params[eBoundTheta]` reads (**242 occurrences across 47 non-test files** — most are chart-agnostic `cov(eBoundPhi, …)` index refs; the ones that matter are the `makeDirectionFromPhiTheta(p[eBoundPhi], p[eBoundTheta])`-shaped reads that would silently misread slopes). Landing this while everything is still spherical is what makes the audit safe: the diff is a no-op for spherical and every reviewed site is already correct before any slope surface exists.

### PR3 — opt-in slope mode on `PlaneSurface`

Add the chart member and constructor argument. Add an integration test for forward tracks through a z-stack with covariance transport enabled. Remove the `1_cm` workaround for slope-mode surfaces only.

### PR4 — `DiscSurface` slope mode

Mirror PR3 on discs. Position stays polar; direction encoding gets the chart switch. Add forward-disc propagation tests.

### PR5 — Telescope example wiring + Python

Flip the example telescope geometry to opt into slope mode. Expose the chart choice in the Python bindings for `PlaneSurface` and `DiscSurface`. Update `TelescopeDetector` example to demonstrate.

### PR6 — Material effects free-space rewrite

Change `computeMaterialEffects` to return `Σ_dir`. Add the sandwich helper. Update `MaterialInteractor`, KF ×3, GSF, and GX2F (deleting the `sin²θ` math). New unit tests for the free-space helper against numerical derivatives on the chart. Regression tests: covariance-transport integration tests should be bit-for-bit invariant for spherical-chart pipelines.

### PR7 — Curvilinear slope-mode

Extend the chart choice to `CurvilinearSurface`. Add `curvilinearChart` to `PropagatorOptions`, default `eSpherical`. Thread it into `boundToCurvilinearTransportJacobian` and the Ridders inline construction sites. Override `CurvilinearSurface::chart()`.

Add covariance-transport tests for forward/backward tracks (`θ ≈ 0, π`) with `eSpherical` (current behavior, marked expected-fail in the regime today) and with `ePlaneSlopes` (should pass cleanly). Remove the `theta close to poles` exclusion in `PropagationEigenConstant.cpp:115` and the `1_cm` shift in `PropagationTests.hpp:306` for the slope-mode runs.

### PR8 — MultiTrajectory accessor helpers and second audit pass

The parameters-object accessors were hardened in PR2.5; this PR extends the same treatment to `MultiTrajectory` proxies (add `state.direction(gctx)` / chart-aware helpers where proxies expose raw bound params) and does a second grep-and-review pass now that slope mode actually exists on planes/discs, catching any callsite that reads raw `eBoundPhi`/`eBoundTheta` off a track state rather than off a parameters object. Must land before PR7 (curvilinear) to bound the curvilinear-audit risk.

### PR9 — IO and EDM4hep boundary

Chart-aware `PrintParameters` names. `RootTrackStatesWriter` branch names and residuals. EDM4hep boundary policy: convert to spherical at write. Documentation update.

### PR10 (optional) — KF/CKF telescope fit benchmark

Demonstrate stable covariance through a barrel-plus-forward-telescope fit, with the barrel running spherical and the telescope arms running slopes (curvilinear states between them in slope mode via `PropagatorOptions`). Use as regression guard.

### Later — flip the curvilinear default

Once PR6, PR7, and PR8 have soaked and the audit is clean, change the default of `PropagatorOptions::curvilinearChart` to `ePlaneSlopes`. This is a behavior change for every pipeline that produces curvilinear states near the z-axis (where the current default silently degrades anyway) and warrants its own release note.

## Open decisions

1. **Accessor semantics.** ~~Keep `phi()/theta()` as "always computed from direction" or deprecate?~~ **Resolved** (see *Accessor semantics*): cache the resolved world direction at construction (resolution A) so accessors stay context-free. Remaining sub-decision: cache 3 doubles vs. recompute-with-stored-transform.
2. **Chart lifetime.** Singletons (recommended — stateless, zero-cost) or per-surface instances (allows charts with parameters)?
3. **Disc inclusion.** Do PR3 and PR4 as separate PRs or one combined "planar surfaces" PR?
4. **Naming.** `PlaneSlopeChart` vs `SlopeChart` vs `CartesianDirectionChart`? Shows up in Python bindings and serialization.
5. **Serialization.** Persist the chart choice in JSON geometry files? If yes, default to `eSpherical` on read for backward compat.
6. **Ordering of PR7 and PR8.** Landing PR8 (accessor audit) first reduces PR7 risk; landing PR7 first delivers the physics fix sooner.

## Risks

- **Context-free `direction()` breaks in slope mode (top risk).** `GenericBoundTrackParameters::direction()`/`phi()`/`theta()` take no `GeometryContext` but slope-mode direction depends on the surface transform. Resolved by caching the world direction at construction (PR2.5, resolution A). If the cache invariant is ever violated — params mutated without re-resolving — every downstream `position()`/`momentum()` silently drifts. Mitigation: make the parameters object own the cache and keep `m_params` mutation funneled through one setter that re-resolves. This is why PR2.5 lands before any slope surface.
- **Index aliasing confusion.** Reusing `eBoundPhi`/`eBoundTheta` for slopes is convenient for the ABI but invites misreads across **242 occurrences in 47 non-test files**. Most are chart-agnostic covariance-index refs; the dangerous shape is `makeDirectionFromPhiTheta(p[eBoundPhi], p[eBoundTheta])`. Mitigations: named accessors that always do the right thing, `eBoundDir0`/`eBoundDir1` aliases for slope-mode code, and the two-pass audit (PR2.5 spherical-only, PR8 after slope exists).
- **Sandwich transpose / wrong chart.** The material projection is `J·Σ·Jᵀ` with `J` the 2×3 free→bound direction Jacobian, and the chart is the *curvilinear* one in the propagator path but the *surface* one in KF/GSF/GX2F. Both are easy to get wrong and neither necessarily fails loudly — a swapped transpose can still compile if `J` is the 3×2 forward Jacobian. Guard: dimension-checked helper signature + the bit-for-bit spherical regression test.
- **Test matrix expands.** Every test that exercises a `PlaneSurface`, `DiscSurface`, or a curvilinear propagation ideally runs in both chart modes. Parametrize the existing tests on the chart.
- **Geometry persistence compatibility.** Old JSON files don't have a chart field; readers must default to `eSpherical`.
- **Downstream callers reading raw `params[eBoundPhi]`.** Existing internal code that does this will silently misread slope-mode parameters. Requires the PR8 audit before slope mode is enabled anywhere production-facing.
- **Curvilinear is everywhere.** PR7's blast radius is wider than the telescope phases because seeders, smoothers, and vertexing all touch curvilinear states. The accessor-recompute rule from PR8 mitigates most of it, but the explicit audit of curvilinear callsites is mandatory before flipping the default later.
- **Material-effects rewrite regression risk.** PR6 changes physics-critical math. Guard: existing covariance-transport integration tests must be bit-for-bit invariant for spherical-chart pipelines. Adds a numerical-derivative unit test for the free-space helper.
- **GSF merging is the trickiest fitter site.** Per-index dispatch on `chart.isPeriodic` keeps the algorithm structure intact but the diff is nontrivial. Isolate in its own commit within PR2 for review.
- **EDM4hep round-trip is lossy.** Slope-mode states convert through spherical at the EDM4hep boundary; the covariance is transported via the chart Jacobian and information is exactly preserved for the direction, but the *representation* changes. Downstream tools reading EDM4hep still see spherical. Document explicitly.

## Out of scope

- **Cylinder slope parametrization.** Doable in principle (slopes in the local tangent frame at the intersection point) but the bound→free Jacobian picks up cross-terms from the position-dependent local frame, and the motivation is weak: tracks that would hit the θ=0 pole geometrically miss a cylinder barrel. Cylinders keep spherical angles.
- **Replacing perigee globally.** Beam-line perigee representation and the vertexing pipeline that consumes it stay spherical. Curvilinear is in scope (PR7); perigee is a separate, larger conversation.
- **Mixed charts on a single surface.** The chart is set per surface (or per propagator option, for curvilinear) at construction and immutable for the lifetime of the surface.
- **Extending the EDM4hep schema.** Not our call; the boundary policy is convert-to-spherical.

## Summary

The change decomposes into a chart policy object hung off `Surface`, a set of small surgical dispatches at the ~10 known angle-aware sites, a free-space rewrite of scattering (which simplifies GX2F as a bonus), and a `PropagatorOptions` knob for curvilinear. The base `Surface` grows exactly one virtual. MultiTrajectory, TrackContainer, KF/CKF/GX2F/GSF *structural* code, and all non-planar surfaces are unchanged.

Two things are **not** free and must not be treated as afterthoughts:

1. **The parameters object stops being self-describing.** `direction()`/`phi()`/`theta()` currently need no `GeometryContext`; slope mode makes the world direction context-dependent. The fix (cache the resolved direction at construction, PR2.5) has to land *before* any slope surface exists, and its cache invariant is the single most correctness-critical piece of the change.
2. **The `eBoundPhi`/`eBoundTheta` reuse is a real audit, not a formality** — 242 occurrences across 47 non-test files, triaged in two passes (spherical-only in PR2.5, post-slope in PR8). Steppers also change: their state must carry the curvilinear chart so `transportCovarianceToCurvilinear` projects scattering correctly.

With those handled, mixed geometries (barrel + telescope) work end-to-end in a single binary.
