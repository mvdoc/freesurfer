# `mris_flatten` regression: FreeSurfer v6 → v7/v8

**Question investigated:** With default parameters, a cortical patch that flattens correctly
under FreeSurfer 6 produces a distorted / failed flat map under FreeSurfer 7 and 8. What
changed, and what is the most likely cause?

**Method:** The local checkout is a shallow clone of the `dev` branch, so the full history
was traced against the upstream tags available in the fork (`v6.0.0` … `v8.2.0`) using the
GitHub API, and cross-checked against the current source tree in this repo.

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

## Hypotheses, ranked

| # | Hypothesis | Confidence | Why |
|---|-----------|-----------|-----|
| **H1** | The `MRISsampleDistances` rewrite + `VERTEX_TOPOLOGY` split changed the sampled `dist_orig` target distances, so the optimizer now minimizes toward a different metric. | **High** | Flattening is *defined* by `dist_orig`; the commit explicitly "fixed bugs" and changed `nsize` semantics; current code uses the `_new` implementation unconditionally. |
| **H2** | The neighborhood-init change in `mris_flatten.c` (`nsize` → `nsizeMax`, plus `...AndDist`) feeds a different neighborhood into distance sampling. | Medium | Directly post-v6, directly in the tool, changes what gets sampled. |
| **H3** | The flatten path now runs with the post-refactor `avg_nbrs`, mis-scaling the spring/smoothing forces — and unlike register/sphere there is no compatibility gate. | Medium | `avg_nbrs` scales the deformation/smoothing terms; the devs flagged it "incorrect" elsewhere but never pinned it for flattening. |
| **H4** | Secondary effects of the refactor: `dist`/`dist_orig` (re)allocation/clear timing, `float`↔`double`, or OpenMP-induced ordering differences. | Low–Medium | "managed xyz/dist" commits touched exactly this; would add noise/instability rather than a clean offset. |

All four share the same root cause: **the 2018–2019 `mrisurf` rewrite between v6.0.0 and
v7.0.0**, not any change to `mris_flatten` arguments or defaults.

---

## Suggested next steps to confirm

1. **Localize the break to a release.** Build `v7.0.0` and rerun the same patch with
   `mris_flatten -norand <in> <out>`. If v7.0.0 already fails, the cause is the 2018 refactor
   (consistent with H1–H3); this rules out anything added in v7.x/v8.x.
2. **Diff the preserved metric directly.** Run with `FS_MEASURE_DISTANCES=1` (handled in
   `MRISunfold`, `mrisurf_integrate.cpp:2396`) under both v6 and v8 on the identical patch and
   diff the resulting `distance.log` (`d` vs `dist_orig`). A systematic difference there
   confirms H1/H2 (the target itself moved); a matching target but different output points to
   H3/H4 (the optimizer/force scaling).
3. **Compare the `flatten.log` distance-error traces** ("starting/final distance error %%")
   between versions to see whether v8 starts from a worse metric or merely converges worse.
4. **Probe `avg_nbrs`.** Add a debug print of `mris->avg_nbrs` just before the `MRISunfold`
   call in both versions; if they differ, H3 is in play. (The `MRIS_SPHERE_NEW_BEHAVIOR` /
   `MRIS_REGISTER_NEW_BEHAVIOR` env vars do **not** affect the flatten path, so they won't fix
   it — but they confirm the developers' awareness of the `avg_nbrs` problem.)
5. **Targeted bisect** of the 2018 refactor PRs in upstream order:
   `92c2723b`/`1b9de0a` (topology split) → `e25a048`/`156a430` (MRISsampleDistances rewrite) →
   `624ffc0a`/`dd421275` (managed xyz/dist) → `9bb0608` (neighborhood reset in the tool).
   H1 predicts the break lands at the `MRISsampleDistances` rewrite.

---

## Key source references (current tree)

- `mris_flatten/mris_flatten.cpp:278-294` — default parameters (unchanged from v6).
- `mris_flatten/mris_flatten.cpp:350-355,387-388` — post-v6 neighborhood/rip changes.
- `mris_flatten/mris_flatten.cpp:511-522` — `MRISreadOriginalProperties` → `MRISflattenPatch` →
  `MRISscaleBrain` → `MRISunfold`.
- `utils/mrisurf_vals.cpp:3484-3689` — rewritten `MRISsampleDistances`/`MRISsampleDistances_new`.
- `utils/mrisurf_integrate.cpp:2283-2587` — `MRISunfold`; `:2376-2384` samples `dist_orig`;
  `:2543-2545` `avg_nbrs` hack.
- `utils/mrisurf_compute_dxyz.cpp:2566,2740` — spring/smoothing terms scaled by `1/avg_nbrs`.
- `utils/mrisurf_metricProperties.cpp:4269` — `MRISsetNeighborhoodSizeAndDist`.
