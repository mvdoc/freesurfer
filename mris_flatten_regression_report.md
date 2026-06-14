# `mris_flatten` regression: FreeSurfer v6 → v7/v8

**Question investigated:** With default parameters, a cortical patch that flattens correctly
under FreeSurfer 6 produces a distorted / failed flat map under FreeSurfer 7 and 8. What
changed, and what is the most likely cause?

**Method:** The local checkout is a shallow clone of the `dev` branch, so the full history
was traced against the upstream tags available in the fork (`v6.0.0` … `v8.2.0`) using the
GitHub API, and cross-checked against the current source tree in this repo.

> **Update (2026‑06‑14): empirical results are in — see the new section below.**
> A side-by-side run of FS6 vs FS8 on the same patch **confirms the regression** but
> **refutes the original top hypothesis (H1)**: the sampled target distances (`dist_orig`)
> are essentially identical between versions, so the cause is on the *optimizer* side, not
> the *target* side. The ranking and the "next steps" below have been revised accordingly,
> and one factual error in the original draft has been corrected (the
> `MRIS_SPHERE_NEW_BEHAVIOR` env gate **does** affect the flatten path).

---

## TL;DR

1. **The `mris_flatten` command-line defaults did not change.** Every default parameter in
   `main()` (`dt`, `tol`, `n_averages`, `l_dist`, `l_nlarea`, `niterations`, `nbhd_size`,
   `max_nbrs`, `scale=3`, `nbrs=2`, `momentum=0.9`, …) is byte-for-byte identical between
   `v6.0.0/mris_flatten/mris_flatten.c` and the current `mris_flatten/mris_flatten.cpp`.
   **The regression is therefore not a parameter change in the tool itself** — it lives in
   the shared surface library (`utils/mrisurf*`).

2. The flattening math was rewritten wholesale between v6 (Jan 2017) and v7 (2020) as part of
   the large **`mrisurf.c` refactor** (Bevin Brett / Andrew Hoopes, 2018–2019): the monolithic
   `mrisurf.c` was split into `mrisurf_*.cpp`, the `VERTEX` struct was split into
   `VERTEX` + `VERTEX_TOPOLOGY`, and the geodesic-distance sampler `MRISsampleDistances` was
   reimplemented.

3. **Primary hypothesis:** the rewrite of `MRISsampleDistances` (and the surrounding
   neighborhood/`dist_orig` bookkeeping) changed the *target geodesic distances* that the
   flattening optimizer tries to preserve. Because flattening is a metric-distortion
   minimization driven entirely by those sampled distances, even a small change in which
   neighbors are sampled or how edge-length distances are accumulated produces a visibly
   different (and sometimes failed) flat map — exactly the reported symptom.

---

## Empirical validation (follow-up run, 2026‑06‑14)

An independent run put the identical patch (sub‑190 rh, continuity-only projection) through
both binaries with bare default parameters, `OMP_NUM_THREADS=8`:
FS6 (`freesurfer-6.0`, Jan 2017) vs FS8 (`freesurfer-8.0.0`, build 8.0.0‑20250204). Both
printed identical defaults, re-confirming the "defaults unchanged" claim.

**Regression confirmed** (independent scoring of the final maps):

| | flipped triangles | Fischl J_d | rel. distortion | opt scale |
|---|---|---|---|---|
| FS6 | 5 / 322,629 (0.002%) | 0.863 mm | 22.2% | 0.915 |
| FS8 | 6,920 / 322,629 (2.145%) | 1.915 mm | 48.8% | 0.587 |

FS8 ends ~1,400× more flipped, ~2.2× the distortion, and over-folded (opt scale 0.59 —
the boundary curls into petals). Mid-run FS8 diverged to **38.7% flipped** before partially
recovering, and took 860 vs FS6's 600 iterations. On a harder patch the 38% excursion may not
recover at all — that is the total-failure mode the user reports.

**H1 is refuted (the key correction).** Using `FS_MEASURE_DISTANCES=1` (which dumps
`distance.log` as `(current_d, dist_orig)` pairs and then exits), over all 16.4 M sampled
pairs the targets (`dist_orig`) are essentially identical between versions:
mean 2.7557 (FS6) vs 2.7553 (FS8); p50/p90/p99 ratios all 0.9998–1.0000. A
`1.207 → 1.09` `TRIANGLE_DISTANCE_CORRECTION` change would have shifted every target ~10%;
the observed global shift is **< 0.02%**. The rewrite does change *which* pairs are sampled
(FS8 has +11,678 pairs, +0.07%) and perturbs ~28% of individual targets by a zero-mean ±1–3%,
but it does **not** move target magnitudes systematically. **Identical targets ⇒ the cause is
the optimizer (H3/H4), not the target (H1/H2).** H2's sampling change is real but far too
small to explain a 2% flip rate.

**Discriminating signature.** In the SSE breakdown, FS8's **area / negative-area energy term
sits persistently ~3× lower than FS6's** (~15.5 vs ~40–50) across the entire run, i.e. FS8
under-penalizes negative area and cannot prevent folding. (Both versions show large
start-of-pass SSE spikes that then anneal — that is normal multi-pass behavior, *not* the
regression; the regression shows up only in the final map quality and the persistent
area-term magnitude gap.) H4 (alloc/clear timing, float/double, OMP ordering) would manifest
as run-to-run noise, not a clean persistent ~3× offset, so it is unlikely to be primary.

### What the code says about the area-energy gap (refinement of H3)

The empirical signature is an **area-term** gap, but the two `1/avg_nbrs` sites
(`mrisurf_compute_dxyz.cpp:2566`, `:2740`) scale the **distance** term, not the area term —
and on the flatten path the relevant `avg_nbrs` looks the same in both versions:

* `MRISsampleDistances_new` **recomputes** `avg_nbrs = total_nbrs / MRISvalidVertices` over the
  *expanded* sampled set (`mrisurf_vals.cpp:4074`), exactly as v6 did — and the measured
  neighbor-count difference is only +0.07%. So the distance-term normalization during the main
  integration epochs is, to first order, **not** different between versions.

Two consequences:

1. The ~3× **area**-energy gap is therefore more plausibly a **symptom** of over-folding than
   its direct cause — the optimizer lets the patch fold, and the (lower) area energy is what a
   folded-but-cheap configuration looks like to the SSE.
2. The one `avg_nbrs` manipulation that *is* reachable on the flatten path **and** is
   version-gated is the final fold-removal block inside `MRISunfold`
   (`mrisurf_integrate.cpp:2543-2547`):
   ```c
   float incorrect_avg_nbrs = mris->avg_nbrs;        // expanded value from sampling
   MRISresetNeighborhoodSize(mris, 1);                // recomputes avg_nbrs to ~1-ring (~6)
   if (getenv("MRIS_SPHERE_NEW_BEHAVIOR") == nullptr) // DEFAULT: restore the old/large value
     mris->avg_nbrs = incorrect_avg_nbrs;
   mrisRemoveNegativeArea(mris, parms, base_averages > 32 ? 32 : base_averages, MAX_NEG_AREA_PCT, 2);
   ```
   This block **runs for flattening** (a flatten is `MRIS_PLANE`, and this is in the common
   body of `MRISunfold`, before the plane-only smoothing at `:2549`). It tries to *emulate*
   the old/v6 behavior by default; if that emulation is imperfect — e.g. `avg_nbrs` entering
   the block differs from what v6 actually had there — this is exactly where a clean,
   persistent force-scaling offset in the fold-removal phase would come from.

**Correction to the original draft:** I previously wrote that `MRIS_SPHERE_NEW_BEHAVIOR` /
`MRIS_REGISTER_NEW_BEHAVIOR` "do not affect the flatten path." That is wrong for
`MRIS_SPHERE_NEW_BEHAVIOR`: it gates lines `2543-2545`, which are inside `MRISunfold` and run
on every flatten. This makes it a **rebuild-free experiment** (see step 1 below).

---

## Background: how `mris_flatten` works

`main()` (`mris_flatten/mris_flatten.cpp`) does, in order:

1. read the surface and the patch/label, rip the rest of the surface;
2. set the neighborhood size and read the *original* (folded) surface properties
   (`MRISreadOriginalProperties`) — this defines the geometry to be preserved;
3. `MRISflattenPatch()` — project the patch onto its average tangent plane (an initial guess);
4. `MRISscaleBrain(...scale=3...)` then `MRISunfold()` — the real work: iteratively move
   vertices in 2-D to minimize the difference between the current inter-vertex distances and
   the **sampled original geodesic distances** (`dist_orig`), plus area/spring/angle terms.

The quantity that defines "correct" is `dist_orig`, produced by **`MRISsampleDistances`**
(called inside `MRISunfold` after restoring the `ORIGINAL_VERTICES` positions —
`utils/mrisurf_integrate.cpp:2376-2384`). If `dist_orig` changes, the optimum changes.

---

## What actually changed (with provenance)

### A. `mris_flatten.c` itself — neighborhood initialization (post-v6)

Commit `9bb0608` *"added resetting of neighborhood size"* (Bruce Fischl, 2018‑11‑21, after
v6.0.0) changed the patch-setup block:

* **v6.0.0:**
  ```c
  MRISresetNeighborhoodSize(mris, mris->vertices[0].nsize) ;
  ```
* **current (v7/v8):** `mris_flatten.cpp:387-388`
  ```c
  MRISsetNeighborhoodSizeAndDist(mris, mris->vertices_topology[0].nsizeMax);
  MRISresetNeighborhoodSize(mris, mris->vertices_topology[0].nsizeMax) ; // set back to max
  ```
  and (lines 350-355) a new block that rips vertices with `vnum == 0`.

Net effect: the surface neighborhood is now forced to `nsizeMax` (= 3) **and distances are
recomputed** (`...AndDist`) before the patch is processed, whereas v6 reset to the surface's
current `nsize`. `MRISsetNeighborhoodSizeAndDist` → `...AndOptionallyDist(..., true)` actively
rebuilds the 2-/3-ring neighbor lists and distance arrays
(`utils/mrisurf_metricProperties.cpp:4269`). This feeds different neighborhoods into the
later distance sampling.

### B. The `VERTEX` / `VERTEX_TOPOLOGY` split and "managed xyz/dist" series (v6 → v7)

Relevant commits (2018):
* `1b9de0a`, `92c2723b` — *"vertex_topology"*: split topology (`v`, `vnum`, `vtotal`, `nsizeCur`,
  `nsizeMax`) out of `VERTEX`.
* `624ffc0a` — *"refurbish mrisurf vals … consolidating who sets VERTEX xyz"*.
* `dd421275` — *"more managed xyz and dist … sometimes clear dist when not valid, begin unifying
  neighbourhood and rip removal"*.

These redefined **when the `dist`/`dist_orig` arrays are (re)allocated, cleared, and considered
valid**, and unified neighborhood vs. rip handling. Any divergence here between the folded-surface
pass and the flattening pass shifts the preserved metric.

### C. `MRISsampleDistances` reimplementation (v6 → v7) — *prime suspect*

Commits:
* `e25a048` — *"Clarify what vertex.nsize means and **fix some bugs** and add a more readable
  MRIsampleDistances implementation"* (Bevin Brett, 2018‑09‑27).
* `3d00ca1` — revert of an earlier attempt, then `156a430` — *"Bevin's tidier mrisample
  distances"* re-landed it.

The current tree dispatches straight into the rewrite: `MRISsampleDistances` is a thin wrapper
around `MRISsampleDistances_new` (`utils/mrisurf_vals.cpp:3484-3506`). This function walks rings
of neighbors out to `nbhd_size` (7 by default), applies a `TRIANGLE_DISTANCE_CORRECTION`
(1.09) edge-length correction, and stores the result in `dist_orig`. The commit message itself
states it **changed the meaning of `nsize` and fixed bugs** — i.e. it deliberately changed the
numbers, with no compatibility flag for `mris_flatten`. This is the most direct route from the
refactor to a different flat map.

### D. The `avg_nbrs` "incorrect value" hacks (no flatten compatibility gate)

During the refactor the developers found that the distance/spring/smoothing terms had been
using an **incorrect `avg_nbrs`** value, and they preserved the *old* behavior behind opt-in
env vars:

* `utils/mrisurf_integrate.cpp:804-808` (`MRISregister` path):
  > "this is an ugly hack to restore the old (and incorrect) computation of the distance term,
  > which had previously used an incorrect value for avg_nbrs"
  gated by `MRIS_REGISTER_NEW_BEHAVIOR` (line 934).
* `utils/mrisurf_integrate.cpp:2543-2545` (`MRISquickSphere` path):
  ```c
  float incorrect_avg_nbrs = mris->avg_nbrs;
  MRISresetNeighborhoodSize(mris, 1);
  if (getenv("MRIS_SPHERE_NEW_BEHAVIOR") == nullptr) mris->avg_nbrs = incorrect_avg_nbrs;
  ```

`avg_nbrs` directly scales the spring/smoothing forces (`norm = 1.0f / mris->avg_nbrs` in
`utils/mrisurf_compute_dxyz.cpp:2566` and `:2740`) and the final fold-smoothing pass in
`MRISunfold` (`mrisurf_integrate.cpp:2549-2566`). Crucially, there is **no `MRIS_FLATTEN_*`
gate** — the flattening path (`MRISunfold`) uses whatever `avg_nbrs` the refactored
neighborhood code now produces. Because `avg_nbrs` is recomputed by the rewritten
`MRISsampleDistances`/`MRISresetNeighborhoodSize` (`mrisurf_topology.cpp:987,2572,2658` and
`mrisurf_vals.cpp:4074`), the flatten path may now be running with a *different* `avg_nbrs`
than v6 did — and unlike register/sphere, nobody pinned it.

### E. Recent dev churn (post‑v8.0.0, probably unrelated to the v6→v8 regression)

* `665f8ed` (2024‑12‑20) — added `-seg`/label-to-patch support.
* `c1b6b2c` (2026‑02‑27, D. Greve) — *"mris_flatten.cpp revert to previous revision from
  2/25/2026"*, i.e. a change made on 2/25 was rolled back two days later.

These post-date the v8 releases the user is running, so they are unlikely to explain a v6→v8
difference, but they show the file is being actively edited and that maintainers are aware of
fragility here.

---

## Hypotheses, ranked (revised after empirical testing)

| # | Hypothesis | Status | Evidence |
|---|-----------|--------|----------|
| **H1** | The `MRISsampleDistances` rewrite moved the sampled `dist_orig` *target* distances. | **Refuted** | `distance.log` targets identical to <0.02% over 16.4 M pairs; a 1.207→1.09 correction would have shifted them ~10%. |
| **H2** | Neighborhood-init / sampling change feeds a different neighbor set. | Real but **insufficient** | Only +0.07% more sampled pairs, zero-mean ±1–3% per-pair jitter — far too small to explain a 2% flip rate. |
| **H3** | The flatten optimizer under-penalizes negative area (force/energy scaling), letting the patch over-fold. | **Leading** | FS8's area/neg-area SSE term persistently ~3× lower; final map is over-folded (opt scale 0.59, petals). Mechanism still to be pinned: the distance-term `1/avg_nbrs` is ~unchanged during the main epochs (`vals.cpp:4074`), so the prime suspect is the **gated fold-removal block** `mrisurf_integrate.cpp:2543-2547` and/or `mrisRemoveNegativeArea` itself. |
| **H4** | Refactor side-effects: alloc/clear timing, `float`↔`double`, OMP ordering. | Unlikely primary | Would produce run-to-run noise, not a clean persistent ~3× offset. |

Net: the regression is **real and reproducible**, lives on the **optimizer** side (negative-area
control during unfolding), and originates in the **2018–2019 `mrisurf` rewrite** between v6.0.0
and v7.0.0 — not in `mris_flatten`'s arguments, defaults, or the sampled target distances.

---

## Suggested next steps to confirm (revised)

Now that H1 is ruled out and H3 leads, the cheapest discriminators come first:

1. **Rebuild-free env-gate test (do this first).** On FS8 run
   `MRIS_SPHERE_NEW_BEHAVIOR=1 mris_flatten -norand <in> <out>` and re-score. This flips the
   `avg_nbrs` used by the fold-removal block (`mrisurf_integrate.cpp:2543-2545`) from the
   "old/large" emulated value to the recomputed 1-ring value. If the flip rate / petal-folding
   changes materially, the regression is localized to that block's `avg_nbrs` handling
   (i.e. the v6 emulation is imperfect). If nothing changes, the negative-area control problem
   is upstream of it — look at `mrisRemoveNegativeArea` and the area-term SSE/gradient
   (`mrisComputeNonlinearAreaSSE` / `mrisurf_sseTerms.cpp`).
2. **Instrument the per-term energies and `avg_nbrs`.** In a source build, print `mris->avg_nbrs`
   plus each SSE component (dist, area, neg-area, spring, angle) at the top of every
   `mrisIntegrationEpoch` in both v6 and v8 on the same patch, and diff. The empirical run
   already shows the ~3× **area** gap; this pins *which* term/coefficient and *when* the two
   versions diverge (main epochs vs the `:2543-2547` fold-removal vs the `:2549` plane smoothing).
3. **Bisect the negative-area path, not the sampler.** With H1 dead, bisect the 2018 refactor in
   upstream order but watch the **area/neg-area** machinery rather than `MRISsampleDistances`:
   `92c2723b`/`1b9de0a` (topology split) → `624ffc0a`/`dd421275` (managed xyz/dist, rip/area
   handling) → `9bb0608` (neighborhood reset in the tool). The `e25a048`/`156a430`
   `MRISsampleDistances` rewrite can be deprioritized (targets proven identical).
4. **Candidate fix to prototype** (once step 1/2 localizes it): add a `MRIS_FLATTEN_*` style
   compatibility gate — or, better, make the fold-removal `avg_nbrs` consistent with the
   neighborhood actually used by `mrisRemoveNegativeArea` — so FS8's negative-area penalty
   matches FS6's during unfolding. Validate by re-scoring flip-rate / J_d / opt-scale against
   the FS6 numbers in the table above.

Artifacts from the follow-up run (both flat patches, both `distance.log`s, scoring/plotting
scripts, side-by-side figure) live at `/data2/projects/autoflatten/fs_regression/`.

---

## Key source references (current tree)

- `mris_flatten/mris_flatten.cpp:278-294` — default parameters (unchanged from v6).
- `mris_flatten/mris_flatten.cpp:350-355,387-388` — post-v6 neighborhood/rip changes.
- `mris_flatten/mris_flatten.cpp:511-522` — `MRISreadOriginalProperties` → `MRISflattenPatch` →
  `MRISscaleBrain` → `MRISunfold`.
- `utils/mrisurf_vals.cpp:3484-3689` — rewritten `MRISsampleDistances`/`MRISsampleDistances_new`.
- `utils/mrisurf_vals.cpp:4074` — `MRISsampleDistances_new` recomputes `avg_nbrs` over the
  expanded sampled set (same as v6) — why the distance-term scaling is ~unchanged in the main epochs.
- `utils/mrisurf_integrate.cpp:2283-2587` — `MRISunfold`; `:2376-2384` samples `dist_orig`;
  `:2543-2547` the `MRIS_SPHERE_NEW_BEHAVIOR`-gated `avg_nbrs` fold-removal hack **that runs on
  the flatten path** (leading suspect for H3).
- `utils/mrisurf_compute_dxyz.cpp:2566,2740` — spring/smoothing terms scaled by `1/avg_nbrs`.
- `utils/mrisurf_metricProperties.cpp:4269` — `MRISsetNeighborhoodSizeAndDist`.
