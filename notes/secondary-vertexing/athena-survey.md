# Secondary vertex finding in Athena — survey, and what it implies for ACTS

**Status:** survey / pre-design. No code changes proposed yet.
**Athena source surveyed:** `/scratch/pagessin/dev/acts-in-athena/athena` (nightly checkout, ~1.9 GB).
**ACTS reference:** this repository, `Core/{include/Acts,src}/Vertexing`, `Examples/Algorithms/Vertexing`.

All Athena paths below are relative to the Athena source root. All ACTS paths are relative to this repository.

---

## 1. Executive summary

Athena has **six** distinct secondary-vertex (SV) reconstruction chains, in four structurally
different families:

| Family | Members | Seeding principle | Output topology |
|---|---|---|---|
| **Pairwise + graph clustering** | `NewVrtSecInclusiveTool`, `VrtSecInclusive`, `InDetVKalVxInJetTool` | fit *every* track pair, build a compatibility graph, extract cliques | N-track, unconstrained |
| **Iterative seed + adaptive fit** | `InDetAdaptiveMultiSecVtxFinderTool` | 3D mode of pairwise DCA midpoints, then AMVF | N-track, unconstrained |
| **Decay chain on an axis** | `InDetSecVxFinderTool` (JetFitter) | tracks assigned to vertices constrained to lie on the jet axis | **multiple vertices on a common flight line** |
| **Exclusive two-track** | `InDetV0Finder`, `InDetConversionFinderTools` | charge-signed pairing + analytic circle intersection | single 2-track vertex, mass-hypothesis tagged |

The dominant, most-developed pattern by far is **pairwise + graph clustering** (Kostyukhin's VKal
family). Three of the six chains are that pattern with different track pre-selection and
post-selection; the algorithmic core is essentially identical between them.

**Provenance** (from Athena git history — first commit adding each package):

| Package | Added | Commits since 2023 | Last touched |
|---|---|---|---|
| `InDetV0Finder` | 2013-12-29 | 19 | 2026-07-13 |
| `VrtSecInclusive` | 2014-04-22 | 41 | 2026-07-23 |
| `InDetSecVxFinderTool` (JetFitter) | 2014-07-19 | 9 | 2026-07-15 |
| `InDetVKalVxInJetTool` | 2014-08-13 | 31 | 2026-07-16 |
| `InDetConversionFinderTools` | 2014-09-03 | 14 | 2026-07-14 |
| **`NewVrtSecInclusiveTool`** | **2020-05-29** | 35 | 2026-06-15 |
| `InDetAdaptiveMultiSecVtxFinderTool` | 2022-07-19 | 11 | 2026-07-21 |

Caveat: the 2013–2014 cluster is the CVS/SVN → git migration, so those dates are an upper bound on
age, not a creation date. The two later dates are genuine.

Two conclusions. **Within the VKal/clique family, `NewVrtSecInclusiveTool` is the newest** — it is
Kostyukhin's rewrite of the same idea, and the reason it reads as the cleanest of the three.
`InDetAdaptiveMultiSecVtxFinderTool` is newer still in absolute terms but belongs to a separate
lineage (its own header dates the work to 2019). And **all six are actively maintained** — none of
this is dead code, so "newest" here means "most recent design", not "the one that replaced the
others".

**The single most important finding for ACTS:** ACTS has a complete, well-tested primary-vertexing
stack (fitters, updaters, impact-point estimation, annealing, finders) but **every one of its seed
finders is one-dimensional in `z`, with `x`/`y` taken from the beam constraint.** See
`Core/src/Vertexing/ZScanVertexFinder.cpp:104-107`. That is the structural blocker: a secondary
vertex is by definition displaced in `x`/`y`, so no existing ACTS seeder can propose it. Everything
else in the ACTS vertexing stack — `AdaptiveMultiVertexFitter`, `FullBilloirVertexFitter`,
`ImpactPointEstimator`, `KalmanVertexUpdater`, `IVertexFinder` — is geometry-agnostic and reusable
as-is.

Consequently the minimum viable baseline is roughly: **one 3D seeder + one selection/merging layer**,
both bolted onto the existing `IVertexFinder` interface. The fitting machinery needs no changes.

---

## 2. Part A — the Athena landscape

### A.1 `NewVrtSecInclusiveTool` (the reference inclusive SV finder)

`Reconstruction/VKalVrt/NewVrtSecInclusiveTool/` (~3.1 kLOC)

The newest member of the family (2020) and the best single template — Kostyukhin's rewrite of the
`VrtSecInclusive` idea, roughly 2.4× smaller than its predecessor for the same core algorithm. Entry
point `findAllVertices()` → `getVrtSecMulti()` in `src/VrtSecMulti.cxx`.

**Pipeline** (`src/VrtSecMulti.cxx:29-529`):

1. **Track selection** — `selGoodTrkParticle()`. Hit-count, χ², pT, |d0|, |z0| cuts; plus an
   anti-pileup cut on the 2D beam–perigee significance.
2. **Two-track vertexing over all pairs** — `select2TrVrt()` in `src/Sel2TrkVertices.cxx:31-239`:
   - per-track 3D IP significance w.r.t. the PV, and a `dR/dZ` significance-ratio cut to kill pileup
     (`src/Sel2TrkVertices.cxx:99-101`);
   - for each surviving pair: **fast crude vertex estimate** (`VKalVrtFitFast`, two-circle
     intersection — see A.7) which also returns the track–track Δz; pairs with large Δz are dropped
     *before* the expensive fit (`:117-118`);
   - full 2-track fit with pion mass hypotheses, χ² < 20 for 1 d.o.f.;
   - a BDT selector (`TwoTrackVrtBDTSelector`) does the real discrimination;
   - optional rejection of vertices sitting on material layers (`distToMatLayerSignificance`);
   - accepted pairs become **edges** in a Boost `adjacency_list` compatibility graph, and the
     crude vertex position is cached per pair (`:213-215`).
3. **Clique extraction** — `bron_kerbosch_all_cliques` over the compatibility graph
   (`src/VrtSecMulti.cxx:90`). Each maximal clique is a candidate N-track vertex; the starting point
   is the *median* of the cached pairwise positions (`estimVrtPos`).
4. **Overlap resolution** — iterative loop (`:140-190`) over vertex pairs sharing tracks, ordered by
   number of shared tracks:
   - if one vertex's track set is contained in another's, kill the subset;
   - if both have fit probability > 1 %, try `mergeAndRefitVertices()`; keep the merge if the merged
     probability > 1 %;
   - otherwise `refineVerticesWithCommonTracks()` — decide which of the two is "better" (by
     multiplicity, then probability) and strip the shared tracks from the worse one.
5. **Proximity merging** — `minVrtVrtDist()` computes normalised vertex–vertex distances; while the
   minimum is below `VertexMergeCut` (default 4σ), merge and refit; on failure drop the
   worst-χ² track and retry (`:230-248`).
6. **Quality repair** — `improveVertexChi2()` iteratively removes the worst-χ² track until the fit
   probability exceeds 1 %; `mostHeavyTrk()` removes the track contributing most to the invariant
   mass if the vertex exceeds `VrtMassLimit`.
7. **Final selection** (`:278-318`):
   - drop vertices *behind* the PV (`projectedVrt < 0`, i.e. momentum–displacement projection);
   - fit probability > `GlobVrtProbCut`;
   - 3D PV–SV distance significance > `SelVrtSigCut`;
   - radius < `MaxSVRadiusCut`;
   - final (tighter) BDT pass on the remaining 2-track vertices.
8. **Output** — full refit for the covariance, decorated with BDT weight, nTracks, charge.

**Default cuts** (`NewVrtSecInclusiveTool/NewVrtSecInclusiveTool.h:149-180`) worth carrying over as
starting points:

| Property | Default | Meaning |
|---|---|---|
| `CutPt` | 500 MeV | track pT |
| `CutD0Max` / `CutD0Min` | 10 / 0 mm | track d0 window |
| `MaxZVrt` | 15 mm | track z impact |
| `TrkSigCut` | 2.0 | track 3D IP significance w.r.t. PV |
| `dRdZRatioCut` | 0.25 | pileup rejection |
| `FastZSVCut` | 10 mm | max track–track Δz from the crude estimate |
| `GlobVrtProbCut` | 0.005 | final vertex fit probability |
| `SelVrtSigCut` | 3.0 | 3D PV–SV distance significance |
| `VertexMergeCut` | 4.0 | vertex–vertex merge significance |
| `MaxSVRadiusCut` | 140 mm | ≈ Pixel outer radius |
| `VrtMassLimit` | 5500 MeV | above this, strip the heaviest track |

There is also a **one-track vertex** mode (`MultiWithOneTrkVrt`, default true): a track detached from
a failing multi-track vertex can be kept as a 1-track "vertex" if the original vertex had > 10 %
probability, subject to double-pT and "crosses ≥ 2 other tracks" requirements
(`src/VrtSecMulti.cxx:199-209`).

### A.2 `VrtSecInclusive` (the older/larger LLP-oriented variant)

`Reconstruction/VKalVrt/VrtSecInclusive/` (~7.5 kLOC)

Same family, more stages, heavier configuration surface. Stage list from
`src/VertexingAlgs.cxx`:

```
extractIncompatibleTrackPairs   :39     all-pairs 2-track fit; builds an *incompatibility* list
findNtrackVertices              :395    Trk::PGraph::pgraphm_ minimal graph covering  (not Bron-Kerbosch)
rearrangeTracks                 :758    move tracks between overlapping vertices
reassembleVertices              :985
associateNonSelectedTracks      :1107   attach tracks that failed the base selection
mergeByShuffling                :1317
mergeFinalVertices              :1467
refitAndSelectGoodQualityVertices :1524
```

Two things here are worth stealing conceptually even if the code is not:

- **`associateNonSelectedTracks`** — a second pass that attaches tracks which failed the *primary*
  track selection to already-found vertices. For displaced vertices this matters a lot, because the
  standard selection is tuned for prompt tracks.
- **PV-compatibility cuts** (`src/VertexingAlgs.cxx:308-325`): angle between the PV→SV displacement
  and the vertex momentum, both in the transverse plane (`vPosMomAngT`) and in 3D
  (`vPosMomAng3D`), plus per-track Δφ consistency. These are cheap, powerful fake killers and are
  reproducible in ACTS with no new infrastructure.

Note this package uses a Fortran-derived graph routine (`Trk::PGraph::pgraphm_`) computing a minimal
covering of the *incompatibility* graph — logically the complement of the Bron–Kerbosch
clique enumeration used by the other two. Same result, worse code. If we port, port the clique form.

### A.3 `InDetVKalVxInJetTool` (b-tagging, jet-directed)

`InnerDetector/InDetRecTools/InDetVKalVxInJetTool/` (~5.3 kLOC)

Structurally identical to A.1 (`src/BTagVrtSecMulti.cxx` mirrors `VrtSecMulti.cxx`, same
Bron–Kerbosch include) but:

- takes a **jet direction** as an extra input and selects tracks within the jet cone;
- explicitly identifies and removes **V0 / material-interaction / conversion** tracks into a
  `TrkFromV0` list before clustering (`src/BTagVrtSecMulti.cxx:108-122`);
- has a dedicated per-track NN classifier (`InDetTrkInJetType`) and a material veto map
  (`InDetMaterialVeto`);
- output is a summary vector (vertex mass, vertex/jet energy fraction, n 2-track vertices, …) rather
  than a plain vertex list.

For a first ACTS baseline this is the *variant*, not the target: it needs jets and a detailed
material map. Worth noting as the eventual b-tagging consumer.

### A.4 `InDetAdaptiveMultiSecVtxFinderTool` (the AMVF-based one — closest to ACTS)

`InnerDetector/InDetRecTools/InDetAdaptiveMultiSecVtxFinderTool/` (~780 LOC — by far the smallest)

This is the most directly relevant precedent, because it is **the primary-vertex AMVF re-pointed at
secondary vertices**. Full loop in `src/InDetAdaptiveMultiSecVtxFinderTool.cxx:124-434`:

```
select tracks  (base selector AND a dedicated SV selector)
repeat up to maxVertices (default 25):
    seed = SeedFinder->findSeed(pvX, pvY, perigees)        # IndexedCrossDistancesSeedFinder
    if seed.z() == 0: break                                # seeder signals "nothing left"
    build a loose constraint vertex at the seed (cov = 1e8 * I)
    for every original track:
        d = ImpactPoint3dEstimator->Estimate3dIP(track, seed)
        doe = d / sqrt(sigma_d0^2 + sigma_z0^2)
        if doe < significanceCutSeeding (default 10): attach to candidate
    if fewer than 2 tracks attached: break
    AdaptiveMultiVertexFitter->addVtxTofit(candidate)
    remove fitted tracks (weight > minWeight) from the seed track pool
    if none removed: remove the track closest to the seed, else abort
    accept if ndf > 0 and >= 2 tracks with weight > minWeight
label 2-track vertices as V0 if they pass a Ks/gamma/Lambda mass window test
```

Two details matter for us:

- The seeder is **`Trk::IndexedCrossDistancesSeedFinder`**, given the PV position as `(vx, vy)`
  reference. This is the 3D seeder ACTS lacks.
- The V0 tagging at the end (`V0check`, `:604-664`) is a pure-kinematics test — three mass
  hypotheses (π⁺π⁻ ≈ K⁰ₛ, e⁺e⁻ ≈ γ, pπ ≈ Λ) plus `a0z < 15 mm` and `Rxy < 500 mm`. About 60 lines,
  no infrastructure, trivially portable.

Caveat: the class is marked `ATLAS_NOT_THREAD_SAFE` and its lifetime management is raw-pointer
heavy. The *algorithm* is the useful part, not the implementation.

### A.5 Exclusive two-track chains

**`InDetV0Finder`** (`InnerDetector/InDetRecAlgs/InDetV0Finder/`, ~1.8 kLOC).
`src/InDetV0FinderTool.cxx` splits tracks by charge, loops over (pos × neg) pairs, gets a starting
point from `VertexPointEstimator::getCirclesIntersectionPoint`, does an unconstrained fit, then runs
**mass-constrained refits** for K⁰ₛ (310), Λ (3122), Λ̄ (−3122) via `TrkV0Fitter`, and applies a
`pointAtVertex` collinearity requirement against the PV. Optional BDT.

**`InDetConversionFinderTools`** (~2.8 kLOC). Same skeleton: `TrackPairsSelector` →
`VertexPointEstimator` → fit → `ConversionPostSelector`, plus a single-track conversion path
(`SingleTrackConversionTool`) for asymmetric conversions where one leg is not reconstructed.

Both are *exclusive* finders: they answer "is this pair a K⁰ₛ / a photon conversion", not "where are
the displaced vertices". They are a natural **second** deliverable, not the baseline.

### A.6 `InDetSecVxFinderTool` — JetFitter (decay chain on the jet axis)

`InnerDetector/InDetRecTools/InDetSecVxFinderTool/` (~3.2 kLOC) plus the fitter in
`Tracking/TrkVertexFitter/TrkJetVxFitter/` (~3.3 kLOC).

Topologically distinct from everything above, and the reason the family count is six rather than
five. Instead of finding independent vertices, JetFitter fits **a set of vertices constrained to lie
on a common flight axis** (the jet axis through the PV), which is the natural description of a
b → c → light decay chain. Components: `JetFitterTrackSelectorTool`, `JetFitterTwoTrackVtxFinderTool`,
`JetFitterV0FinderTool`, `JetFitterMultiStageFit`, `InDetImprovedJetFitterVxFinder`; EDM in
`Tracking/TrkEvent/VxJetVertex/`.

It is in production b-tagging (`PhysicsAnalysis/JetTagging/JetTagAlgs/BTagging/python/JetSecVtxFindingAlgConfig.py`,
`BTagLightSecVertexingConfig.py`, `JetFitterSequentialVertexFitterConfig.py`).

**Out of scope for a baseline** — it needs jets, a flight-axis constraint in the fitter, and an
EDM for multi-vertex chains, none of which ACTS has. Recorded here so the taxonomy is complete and
so that "secondary vertexing" is not silently equated with "one displaced vertex at a time".

### A.7 Shared low-level primitives

These are the reusable atoms. Sizes are small; all are portable.

| Primitive | Location | ~LOC | What it does |
|---|---|---|---|
| `vkvFastV` | `Tracking/TrkVertexFitter/TrkVKalVrtCore/src/VKvFast.cxx:42-195` | 150 | Analytic two-circle intersection of two perigees → 3D point; **returns the track–track \|Δz\| as a quality measure**. Handles nested, intersecting and disjoint circles. Used everywhere as the crude 2-track seed. |
| `VertexPointEstimator::intersectionImpl` | `InnerDetector/InDetRecTools/InDetConversionFinderTools/src/VertexPointEstimator.cxx:83+` | 400 | CTVMFT-derived variant of the same idea, with configurable ΔR/ΔZ/R windows and arc-length constraints, per "flag" category. |
| `CrossDistancesSeedFinder` | `Tracking/TrkVertexFitter/TrkVertexSeedFinderTools/src/CrossDistancesSeedFinder.cxx` | 260 | For all track pairs, compute the true 3D DCA midpoint, weight it by `1/(cutoff+d)^p` and by a sigmoid on the constraint χ², then take the **3D mode** of the weighted point cloud. |
| `Mode3dTo1dFinder` | `Tracking/TrkVertexFitter/TrkVertexSeedFinderUtils/src/Mode3dTo1dFinder.cxx` | 110 | 3D mode = independent FSMW 1D mode in x, y, z. |
| `Mode3dFromFsmw1dFinder` | same dir | 1270 | The heavyweight version: clusters in 3D, returns per-seed track index lists via `IMode3dInfo`. Used by `IndexedCrossDistancesSeedFinder`. |
| `NewtonTrkDistanceFinder` / `Trk2dDistanceSeeder` / `SeedNewtonTrkDistanceFinder` | same dir | 290 / 440 / 100 | Exact track–track 3D minimum distance: 2D circle-crossing seed, then Newton iteration on the two arc-length parameters. |

**`Mode3dTo1dFinder` is a ~50-line job in ACTS** because ACTS already has the FSMW 1D mode finder:
`Core/include/Acts/Vertexing/FsmwMode1dFinder.hpp` (`getMode(std::vector<std::pair<double,double>>)`).
It is currently used only inside `ZScanVertexFinder`.

### A.8 EDM and truth matching

- **`Tracking/TrkEvent/VxSecVertex/`** — `VxSecVertexInfo` and derivatives (`VxSecVKalVertexInfo`,
  `VxJetFitterVertexInfo`). The SV finders return a *container of vertices plus metadata*, not a bare
  vertex list. ACTS's `IVertexFinder::find()` returns `std::vector<Vertex>` — likely sufficient for a
  baseline, but note the asymmetry.
- **`InDetSecVtxTruthMatchTool`** — this is the validation contract, and it is worth adopting
  wholesale. Definitions from
  `InnerDetector/InDetRecTools/InDetSecVtxTruthMatchTool/InDetSecVtxTruthMatchTool/InDetSecVtxTruthMatchTool.h`:

  *Reco vertex classes:* `Matched` (> 50 % of track weight from one truth interaction), `Merged`,
  `Split`, `Fake`, `Other`.

  *Truth vertex classes (nested subsets):* `Reconstructable` (fiducial + ≥ 2 charged daughters) ⊃
  `Accepted` (≥ 2 reco tracks) ⊃ `Seeded` (tracks pass cuts) ⊃ `Reconstructed` ⊃
  `ReconstructedSplit`.

  It also carries an origin taxonomy (`KshortDecay`, `LambdaDecay`, `GammaConversion`,
  `HadronicInteraction`, `BHadronDecay`, `DHadronDecay`, `Pileup`, `Fake`, …) which is exactly the
  breakdown we would want in an ACTS performance writer.

---

## 3. Part B — purpose topology

Six chains exist not because Athena is disorganised but because "secondary vertex" is six different
physics objects. This section maps **purpose → algorithm → consumer**, established from the
configuration layer rather than from the source comments.

### B.1 Purpose → algorithm → consumer

| Physics purpose | Algorithm | Configured in |
|---|---|---|
| **b/c-hadron decay vertices in jets** (b-tagging) | `InDetVKalVxInJetTool` | `JetTagAlgs/BTagging/python/JetSecVtxFindingAlgConfig.py`, `BTagLightSecVertexingConfig.py` |
| **b→c decay chain on the jet axis** (b-tagging) | JetFitter (`InDetSecVxFinderTool`) | `JetSecVtxFindingAlgConfig.py`, `JetFitterSequentialVertexFitterConfig.py` |
| **Soft secondary vertices away from jets** (SSV, pileup/b-content weights) | `NewVrtSecInclusiveTool`, WPs `NVSI_SecVrt_{Tight,Medium,Loose}` | `FTagAnalysisAlgorithms/python/SSVWeightsAlgConfig.py` |
| **Long-lived particles / displaced vertices** | `VrtSecInclusive` (×6 tunings) + `NewVrtSecInclusiveTool` | `DerivationFrameworkLLP/python/LLP1.py`, `DerivationFrameworkSUSY/python/SUSY20.py` |
| **Tracking/material studies, V0 & DV in ID** | `NewVrtSecInclusiveTool` (3 presets) + `InDetV0Finder` | `DerivationFrameworkInDet/python/IDTR2.py` |
| **K⁰ₛ / Λ reconstruction** (B-physics) | `InDetV0Finder` + `TrkV0Fitter` | `DerivationFrameworkBPhys/python/V0ToolConfig.py`, `BPHY24.py` |
| **Photon conversions** | `InDetConversionFinderTools` | `egammaAlgs/python/EMVertexBuilderConfig.py` |
| **B-physics trigger** (J/ψ, Bmumux…) | `InDetConversionFinderTools` primitives | `TrigBphysHypo/python/Trig*ComboHypoConfig.py` |
| **LLP trigger** | `VrtSecInclusive` | `TrigTools/TrigVrtSecInclusive/` |
| **ID standalone SV reco** | `InDetAdaptiveMultiSecVtxFinderTool` | `InDetConfig/python/InDetSecVtxFinderConfig.py` — **see B.3** |

### B.2 One algorithm, many purpose presets

This is the most transferable observation in the whole survey. The *algorithm* is generic; the
*purpose* lives entirely in the cut tuning. `NewVrtSecInclusiveTool` ships **seven named presets** in
`Reconstruction/VKalVrt/NewVrtSecInclusiveTool/python/NewVrtSecInclusiveConfig.py`:

| Preset | `CutPt` | `SelVrtSigCut` | `MaxSVRadiusCut` | `VrtMassLimit` | Purpose |
|---|---|---|---|---|---|
| `SoftBFinderTool` | 500 | 2.5 | 50 mm | — | soft b/c vertices |
| `InclusiveBFinderTool` | 500 | 3.0 | 140 mm (default) | — | inclusive b/c |
| `HighPtBFinderTool` | 1000 | 3.0 | 140 mm (default) | — | high-pT b/c |
| `MaterialSVFinderTool` | 500 | 10.0 | 140 mm (default) | 8 000 | material interactions |
| `KsFinderTool` | 1000 | 8.0 | 350 mm | 800 000 | K⁰ₛ |
| `DVFinderTool` | 1000 | 8.0 | 350 mm | 1 000 000 | displaced vertices / LLP |
| `V2TCalibrationTool` | 400 | 2.0 | 140 mm | 5 500 | 2-track BDT calibration |

The pattern is legible: **displacement significance and radius acceptance scale with lifetime**
(2.5σ / 50 mm for soft-b, 8σ / 350 mm for K⁰ₛ and DV), and the mass limit is effectively disabled
for the "find anything displaced" presets.

`VrtSecInclusive` is used the same way — LLP1 instantiates it **six times** in one job
(`DerivationFrameworkLLP/python/LLP1.py:348-431`): nominal, short-lifetime, short-lifetime-no-d0
(for the LRSM dHNL analysis), disappearing-track + large-radius-tracking, plus two track-systematics
variants. IDTR2 instantiates `NewVrtSecInclusiveTool` three times as `_Material`, `_DV`, `_Ks`.

**Implication for ACTS:** the deliverable is not "a secondary vertex finder" but "a finder plus a
presets mechanism". Getting the configuration surface right — which cuts are exposed and in what
units — matters as much as the algorithm, and the seven presets above are a ready-made validation
matrix.

### B.3 Correction: the AMVF-based tool is not in a production chain

Worth stating plainly because it changes the weight of the recommendation in Part E.
`InDetAdaptiveMultiSecVtxFinderTool` is reachable only through `InDetSecVtxFinderAlgCfg`, whose sole
caller in the entire repository is
`InnerDetector/InDetRecAlgs/InDetSecVtxFinder/share/runInDetSecVtxFinder.py` — a standalone runner
script. It appears in no derivation, no reconstruction chain, and no trigger configuration. It is
also marked `ATLAS_NOT_THREAD_SAFE`.

So it is best read as **a well-formed reference implementation of the AMVF-for-SV idea, not a
validated production algorithm**. It remains the closest structural match to ACTS and the cheapest
path to a baseline, but "mirrors a production Athena algorithm" is not an argument that can be made
for it — that argument belongs to the VKal family.

---

## 4. Part C — the pattern, distilled

Stripping framework noise, every Athena inclusive SV finder is the same seven-stage pipeline:

```
1. TRACK SELECTION       quality + displacement (IP significance w.r.t. PV)
2. PAIR SEEDING          all pairs -> analytic crude vertex -> cheap reject -> 2-track fit
3. PAIR SELECTION        chi2 / mass / material / PV-collinearity / (BDT)
4. CLUSTERING            compatibility graph -> maximal cliques -> N-track candidates
5. OVERLAP RESOLUTION    contained-in / merge-and-refit / strip-shared-tracks
6. PROXIMITY MERGING     vertex-vertex significance -> merge -> refit
7. FINAL SELECTION       prob, PV-SV significance, radius, mass, momentum-projection sign
```

The AMVF variant (A.4) replaces stages 2–6 with `3D seed → attach compatible tracks → adaptive fit →
remove used tracks → repeat`, which is exactly the shape of ACTS's existing
`AdaptiveMultiVertexFinder::find()`.

**That is the key insight for the plan:** ACTS's `AdaptiveMultiVertexFinder` already implements the
A.4 control flow. Give it a 3D seeder and appropriate cuts and it becomes a secondary vertex finder.
Stage 1 and stage 7 are then the only genuinely new logic, and they are cut-application, not
algorithms.

---

## 5. Part D — ACTS current state and gap analysis

### D.1 What already exists and is directly reusable

`Core/include/Acts/Vertexing/`:

| Component | Reusable for SV? | Note |
|---|---|---|
| `AdaptiveMultiVertexFitter` / `-Finder` | **yes, unchanged** | finder control flow already matches A.4 |
| `IterativeVertexFinder` | yes, unchanged | alternative host for a 3D seeder |
| `FullBilloirVertexFitter` | yes, unchanged | natural fitter for small N-track candidates |
| `KalmanVertexUpdater` (+ Impl3/Impl4) | yes, unchanged | |
| `ImpactPointEstimator` | yes | has `calculateDistance(track, vtxPos)`, `getVertexCompatibility`, `getImpactParameters`, `estimate3DImpactParameters` — i.e. the ACTS equivalent of `ImpactPoint3dEstimator`, already there |
| `HelicalTrackLinearizer`, `NumericalTrackLinearizer` | yes, unchanged | |
| `AnnealingUtility` | yes, unchanged | |
| `FsmwMode1dFinder` | yes | the 1D building block for a 3D mode finder |
| `Vertex`, `TrackAtVertex`, `InputTrack`, `VertexingOptions`, `IVertexFinder` | yes | EDM/interfaces are geometry-agnostic |
| `ZScan` / `TrackDensity` / `GridDensity` / `AdaptiveGridDensity` VertexFinder | **no** | all 1D-in-z; see D.2 |
| `HoughVertexFinder` | partially | 3D-ish but designed for a single primary vertex in heavy-ion events |

### D.2 The blocker

`Core/src/Vertexing/ZScanVertexFinder.cpp:104-107`:

```cpp
// constraint x()/y() equals 0 if no constraint
Vector4 output(vertexingOptions.constraint.position().x(),
               vertexingOptions.constraint.position().y(), ZResult,
               vertexingOptions.constraint.time());
```

The same holds for the grid/density seeders: they histogram `z0` (and optionally `t`) and return the
beam-line `x`,`y`. Every ACTS vertex seed therefore lies on the beam axis. There is currently **no
way to propose a candidate at `r > 0`**.

### D.3 Gap table

| Needed | Athena reference | ACTS status | Rough size |
|---|---|---|---|
| Analytic 2-track crude vertex (two-circle intersection + Δz) | `VKvFast.cxx` | **missing** | ~200 LOC + tests |
| Exact track–track 3D DCA | `NewtonTrkDistanceFinder`, `Trk2dDistanceSeeder` | **missing** (`ImpactPointEstimator` does track↔*point*, not track↔track) | ~300 LOC + tests |
| 3D mode of a weighted point cloud | `Mode3dTo1dFinder` | **missing**, but `FsmwMode1dFinder` exists | ~60 LOC |
| Cross-distances 3D seed finder (`IVertexFinder`) | `CrossDistancesSeedFinder` | **missing** | ~250 LOC |
| Compatibility-graph clique clustering | Boost `bron_kerbosch_all_cliques` | **missing** (Boost.Graph is not currently an ACTS Core dependency — needs a decision) | ~150 LOC or dependency |
| Vertex merge / overlap resolution utilities | `MultiUtilities.cxx` | partially — AMVF has `isMergedVertex`, but no track-set overlap resolution | ~250 LOC |
| SV quality selection (PV–SV significance, momentum projection sign, radius, mass) | `VrtSecMulti.cxx:278-318` | **missing** | ~150 LOC |
| V0 kinematic tagging | `V0check` | **missing** | ~80 LOC |
| SV truth matching + performance writer | `InDetSecVtxTruthMatchTool` | **missing** | Examples-level |
| Mass-constrained vertex fit | `TrkV0Fitter` | **missing** | out of scope for a baseline |

---

## 6. Part E — candidate baselines

Three options, in increasing order of new code.

### Option 1 — "AMVF + 3D seeder" (recommended baseline)

Add a `CrossDistancesSeedFinder` implementing `Acts::IVertexFinder`, plus the two primitives it needs
(track–track 3D DCA; 3D FSMW mode). Feed it to the **existing** `AdaptiveMultiVertexFinder` with
SV-appropriate configuration (no beam constraint, `useSeedConstraint`, loose
`tracksMaxZinterval`). Add a thin `SecondaryVertexSelector` applying the stage-7 cuts.

- Structurally identical to an existing Athena tool (A.4), so the design is not speculative.
- Touches no existing ACTS code paths; purely additive.
- Reuses the fitter, updater, IP estimator, annealing entirely.
- **Provenance caveat (B.3):** that Athena tool is *not* part of any production chain — it is a
  standalone reference implementation. So this option inherits a sound structure but no physics
  validation record. The tuning has to be established by us, not inherited.
- Biggest risk: `AdaptiveMultiVertexFinder`'s internals assume beamline-ish geometry in places
  (`tracksMaxZinterval`, the merge-significance logic when `doFullSplitting == false` is z-only).
  Needs auditing, and `doFullSplitting = true` is probably mandatory.

### Option 2 — "Pairwise + cliques" (the VKal pattern)

Port stage 2–7 of `NewVrtSecInclusiveTool` on top of `FullBilloirVertexFitter`.

- Closest to what ATLAS actually uses for inclusive/LLP SV, and to what b-tagging expects — and,
  unlike Option 1, it carries a real production validation record (B.1).
- Highest fidelity for low-multiplicity displaced vertices.
- Comes with a ready-made preset matrix (B.2) to validate against.
- More new code, and O(N²) pair fits — needs a cheap pre-filter (that is precisely what `vkvFastV`
  is for).
- Requires a clique enumerator decision (vendor a small Bron–Kerbosch, or take Boost.Graph).

### Option 3 — exclusive two-track (V0 / conversion)

Charge-split pairing + analytic seed + 2-track fit + mass-window tagging.

- Smallest, most testable, cleanest validation story (K⁰ₛ mass peak is a hard, unambiguous metric).
- But it is *not* general secondary vertexing.

**Suggested shape of the work:** build the shared primitives (two-circle seed, track–track DCA, 3D
mode) first — they are needed by all three options — then Option 1 as the baseline, with Option 3 as
the validation vehicle because it gives an unambiguous physics metric early. Option 2 becomes the
follow-up once the primitives are proven.

The B.3 finding weakens the case for Option 1 somewhat but does not overturn it: it is still the
cheapest route to *something that runs*, and the primitives it needs are shared with Option 2
regardless. But if the goal is a baseline that can be defended as "this is what ATLAS does", that is
Option 2, and the extra cost is the clique layer plus the overlap/merging utilities — call it
+400 LOC over Option 1. Worth an explicit decision rather than defaulting to the cheaper one.

---

## 7. Part F — validation surface

What exists in ACTS today: `Examples/Io/Root/src/RootVertexNTupleWriter.cpp`,
`RootVertexWriter.cpp`, `RootVertexReader.cpp`; `Examples/Scripts/Python/vertex_fitting.py`;
`Python/Examples/src/Vertexing.cpp` for bindings;
`Tests/UnitTests/Core/Vertexing/` with per-component tests and `VertexingDataHelper.hpp`.

What a baseline needs:

1. **Unit tests for the new primitives.** The two-circle intersection and the track–track DCA are
   analytically checkable against straight-line and known-helix configurations — high-value,
   low-cost tests.
2. **A truth-matching layer** following the `InDetSecVtxTruthMatchTool` definitions (A.8). The
   `Reconstructable → Accepted → Seeded → Reconstructed` chain is the right efficiency
   denominator ladder; without it, "efficiency" is meaningless for displaced vertices.
3. **A physics metric.** K⁰ₛ → π⁺π⁻ invariant-mass peak position and width from ttbar or minimum-bias
   samples. This is the fastest way to know the thing works at all.
4. **A fake-rate metric** binned in radius — displaced fakes are radius-dependent (material layers),
   and Athena's material veto exists precisely because of that.

Open question: which sample. The Open Data Detector has no realistic material map for
LLP-style studies, but is fine for K⁰ₛ and for geometric efficiency.

---

## 8. Part G — open questions for the planning phase

1. **Physics target.** Inclusive displaced vertices (LLP-style), b-tagging SV in jets, or V0
   reconstruction? These have materially different track selections and cut tunings. The Athena code
   splits along exactly this axis (Part B). *This is the decision that gates everything else.*
   Note that Part B also shows the target may be less binding than it looks: one algorithm with a
   preset mechanism covers soft-b, inclusive-b, material, K⁰ₛ and DV in Athena. Choosing "which
   preset first" is a lighter decision than choosing "which algorithm".
2. **Core or Examples?** A `CrossDistancesSeedFinder` belongs in `Core/Vertexing`. A jet-directed
   finder arguably does not (it needs jets).
3. **Boost.Graph.** Option 2 needs clique enumeration. Is a Boost.Graph dependency in Core
   acceptable, or do we vendor ~150 lines of Bron–Kerbosch?
4. **Does `AdaptiveMultiVertexFinder` survive off-axis use?** Needs an audit of `tracksMaxZinterval`,
   `doFullSplitting`, `maxMergeVertexSignificance` and the seed-constraint covariance handling before
   committing to Option 1.
5. **PV as input.** Every Athena SV finder consumes a primary vertex (for IP significance, for the
   momentum-projection sign, for the seeder reference point). The ACTS algorithm needs the same
   input contract — this affects the `Examples` wiring and the `IVertexFinder` usage.
6. **Multi-track vs two-track scope for v1.** Restricting the baseline to 2-track vertices removes
   stages 4–6 entirely and makes the first version dramatically smaller. Worth considering as a
   deliberate scope cut.
7. **Is a preset mechanism in scope for v1?** B.2 argues the configuration surface is half the
   deliverable. Deciding early whether presets are a v1 feature or a later addition affects how the
   `Config` structs are shaped, and retrofitting that is more painful than designing for it.
8. **Do we care about decay chains at all?** JetFitter (A.6) is a genuinely different output
   topology — several vertices constrained to one flight axis — and nothing in the proposed baseline
   moves toward it. Fine to declare out of scope, but it should be declared rather than overlooked.

---

## Appendix — file index

**Athena** (relative to Athena root):

```
Reconstruction/VKalVrt/NewVrtSecInclusiveTool/          inclusive SV, cleanest reference   (A.1)
    src/VrtSecMulti.cxx                                 main pipeline
    src/Sel2TrkVertices.cxx                             pairwise seeding
    src/MultiUtilities.cxx                              merge/refine/overlap utilities
    NewVrtSecInclusiveTool/NewVrtSecInclusiveTool.h      config defaults
Reconstruction/VKalVrt/VrtSecInclusive/                  older LLP-oriented variant        (A.2)
    src/VertexingAlgs.cxx                               8-stage pipeline
InnerDetector/InDetRecTools/InDetVKalVxInJetTool/        b-tagging, jet-directed           (A.3)
    src/BTagVrtSecMulti.cxx
InnerDetector/InDetRecTools/InDetAdaptiveMultiSecVtxFinderTool/   AMVF-based               (A.4)
    src/InDetAdaptiveMultiSecVtxFinderTool.cxx
InnerDetector/InDetRecAlgs/InDetV0Finder/                V0                                (A.5)
    src/InDetV0FinderTool.cxx
InnerDetector/InDetRecTools/InDetConversionFinderTools/  conversions                       (A.5)
    src/VertexPointEstimator.cxx                        circle-intersection seed
InnerDetector/InDetRecTools/InDetSecVxFinderTool/        JetFitter decay chain             (A.6)
Tracking/TrkVertexFitter/TrkJetVxFitter/                 JetFitter fitter + EDM helpers    (A.6)
Tracking/TrkEvent/VxJetVertex/                           JetFitter EDM                     (A.6)
Tracking/TrkVertexFitter/TrkVKalVrtCore/src/VKvFast.cxx  crude 2-track vertex              (A.7)
Tracking/TrkVertexFitter/TrkVertexSeedFinderTools/src/CrossDistancesSeedFinder.cxx          (A.7)
Tracking/TrkVertexFitter/TrkVertexSeedFinderUtils/src/Mode3dTo1dFinder.cxx                  (A.7)
Tracking/TrkVertexFitter/TrkVertexSeedFinderUtils/src/NewtonTrkDistanceFinder.cxx           (A.7)
Tracking/TrkEvent/VxSecVertex/                           SV EDM                            (A.8)
InnerDetector/InDetRecTools/InDetSecVtxTruthMatchTool/   validation definitions            (A.8)
Tracking/Acts/ActsVertexReconstruction/                  existing ACTS-in-Athena PV tools
```

**ACTS** (relative to this repository):

```
Core/include/Acts/Vertexing/       IVertexFinder, AdaptiveMultiVertexFinder/-Fitter,
                                   IterativeVertexFinder, FullBilloirVertexFitter,
                                   ImpactPointEstimator, KalmanVertexUpdater,
                                   FsmwMode1dFinder, ZScanVertexFinder, Vertex, TrackAtVertex
Core/src/Vertexing/                implementations
Examples/Algorithms/Vertexing/     AdaptiveMultiVertexFinderAlgorithm, IterativeVertexFinderAlgorithm,
                                   VertexFitterAlgorithm, TruthVertexSeeder
Examples/Io/Root/src/              RootVertexNTupleWriter, RootVertexWriter, RootVertexReader
Python/Examples/src/Vertexing.cpp  python bindings
Tests/UnitTests/Core/Vertexing/    unit tests
```
