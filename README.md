# sPHENIX Solenoid Field Map — CERN Mapping 2022-12-02

Analysis of the sPHENIX solenoid magnetic-field measurement from the CERN
mapping campaign (geometry-corrected point cloud, **2022-12-02**), together with
a direct comparison against the calculated **OPERA** field map used in sPHENIX
offline tracking.

This is a single, self-contained repository covering both the measured-map
analysis and the OPERA comparison. Going forward only two repositories are
relevant for the sPHENIX solenoid field:

- this repository — the **measured** map and its analysis, and
- [`haggerty/sphenixoperamaps`](https://github.com/haggerty/sphenixoperamaps) —
  the **calculated** (OPERA) map.

> ### Supersedes `sphenix-cernfinal-map`
> This repository replaces the earlier
> [`haggerty/sphenix-cernfinal-map`](https://github.com/haggerty/sphenix-cernfinal-map)
> and [`haggerty/cern-opera-comparison`](https://github.com/haggerty/cern-opera-comparison).
> Those were built from an earlier map whose solenoid magnetic centre sat at
> z ≈ −240 mm — **~265 mm** from where the corrected 2022-12-02 point cloud places
> it (**z ≈ +26 mm**, in agreement with OPERA's ≈ +28.5 mm design coil offset).
> *How the earlier and 2022-12-02 maps relate — a survey reanalysis vs. a separate
> correction — is to be confirmed with the CERN mapping group.* The old
> repositories should be considered archival.

## Data files

The original 2022-12-02 mapping tarball **`download.tar` is committed at the repo
root** (~48 MB), so the measured-field inputs travel with the analysis. Unpack it
into `data/` (git-ignored) before running:

```bash
mkdir -p data && tar xf download.tar -C data
```

This gives the repo everything needed to reproduce the **measured** map and all
measured-field analysis. The **OPERA** comparison and the `opera_matched_to_mapping`
delivered map additionally need the OPERA map, which is **not** committed (it
belongs to [`sphenixoperamaps`](https://github.com/haggerty/sphenixoperamaps) /
the `FIELDMAP_TRACKING` CDB) — place it at `data/sphenix3dmapxyz.root`
(±80/±100 cm comparison file) and/or `data/sphenix3dtrackingmapxyz.root` (111³
production grid). The two maps' grids and sizes, and how to locate the tracking
map in the CDB (MD5 content hash, CVMFS path, and the global-tag/IOV lookup),
are documented in [the tracking-field-map section of
`sphenixoperamaps`](https://github.com/haggerty/sphenixoperamaps#the-tracking-field-map).
`sphenix3dmapxyz.root` is exactly the central core of the 111³ tracking map.

| File in `download.tar` | Size | Description |
|------------------------|------|-------------|
| `pointCloudFineFullField.csv`  | 35.8 MB | full field, 2 cm step, 10° azimuthal (~200 k points) |
| `pointCloudRoughFullField.csv` |  7.5 MB | full field, 10 cm step, 10° azimuthal (~42 k points) |
| `pointCloudRoughHalfField.csv` |  7.5 MB | **half** field, 10 cm step, 10° azimuthal (not used here) |

On SDCC the tarball and unpacked CSVs are also at
`/sphenix/data/data02/sphenix/MagnetMapping/cern_2022-12-02/`.

| File in `download.tar` | Size | Description |
|------------------------|------|-------------|
| `pointCloudFineFullField.csv`  | 35.8 MB | full field, 2 cm step, 10° azimuthal (~200 k points) |
| `pointCloudRoughFullField.csv` |  7.5 MB | full field, 10 cm step, 10° azimuthal (~42 k points) |
| `pointCloudRoughHalfField.csv` |  7.5 MB | **half** field, 10 cm step, 10° azimuthal (not used here) |

**CSV format** (surveyor frame, mm and Tesla):

```
x_s, y_s, z_s, |B|, Bx_s, By_s, Bz_s
```

**Surveyor → sPHENIX transform** (applied in `sPHENIXFieldMap`):
`x_phx = x_s`, `y_phx = z_s`, `z_phx = −y_s`;
`Bx_phx = Bx_s`, `By_phx = Bz_s`, `Bz_phx = −By_s`.

**Measured extent** of the raw points, in sPHENIX coordinates (the data
coverage, distinct from the interpolation grid the loader builds):

| Map | points | x (mm) | y (mm) | z (mm) | r (mm) | φ |
|-----|-------:|--------|--------|--------|--------|---|
| `pointCloudFineFullField`  | ~200 k | −859 … 851 | −857 … 853 | −1784 … 2233 | 52 … 859 | full 360° |
| `pointCloudRoughFullField` | ~42 k  | −859 … 851 | −857 … 853 | −2192 … 2206 | 52 … 859 | full 360° |

- **Transverse:** cylinder of radius **~860 mm**, full azimuth, with a **central
  hole** — the smallest sampled radius is **~52 mm** (nothing is measured on axis).
- **Axial:** the fine map spans **z ≈ −1.78 … +2.23 m**; the rough map reaches
  ~400 mm further south (**z down to ≈ −2.19 m**). The rough map's only unique
  contribution is this −z (south) fringe tail below ≈ −1785 mm — elsewhere the
  fine map overrides it.

The OPERA map is a `TNtuple` named `fieldmap` with branches `x,y,z` (cm) and
`bx,by,bz` (T) on a uniform Cartesian grid.

## Repository layout

```
.
├── sPHENIXFieldMap.{h,cxx}   # field-map class: reads CSVs, phi-averages onto an
│                             #   (r,z) grid, enforces ∇·B = 0 for a consistent Br
├── analysis/
│   ├── checkFieldMap.C       # Bz(r,z), Br(r,z), ∇·B(r,z), Bz on axis; Maxwell check
│   ├── findCenter.C          # magnetic centre (Bz peak / Br sign change)
│   ├── checkAlignment.C      # solenoid-axis tilt vs sPHENIX z (global φ-averaged m=0)
│   └── estimateTilt.C        # two independent tilt estimators (near-axis + axis-line fit)
├── comparison/
│   ├── OperaMap.h            # shared OPERA loader/interpolator (grid read from ntuple)
│   ├── compareNewVsOpera.C   # quick overlays measured vs OPERA
│   ├── compareFieldMaps.C    # full m=0/m=1 Fourier comparison + Maxwell residuals
│   ├── findCenterOpera.C     # OPERA magnetic centre (same estimators as findCenter.C)
│   └── checkTilt.C           # azimuthal m=1 tilt signature from the raw CSVs
├── export/
│   └── makeDeliveredMaps.C  # the two PHField3DCartesian drop-in maps for reco
├── plots/                    # generated PDFs/PNGs (committed)
└── data/, output/            # git-ignored (inputs / large derived ROOT files)
```

All macros are run from the repository root and compiled with ACLiC (`+`); each
takes a data directory and/or output directory argument.

## How to run

```bash
# Field-map analysis                                      # plots written to plots/
root -l -b -q 'analysis/checkFieldMap.C+("data","plots")'   # -> fieldMap_overview.{pdf,png}
root -l -b -q 'analysis/findCenter.C+("data")'              # console only (centre estimators)
root -l -b -q 'analysis/checkAlignment.C+("data","plots")'  # -> fieldMap_alignment.{pdf,png}
root -l -b -q 'analysis/estimateTilt.C+("data","plots")'    # -> tilt_estimates.{pdf,png}

# OPERA comparison
root -l -b -q 'comparison/compareNewVsOpera.C+("data","data/sphenix3dmapxyz.root","plots")'  # -> compare_newVsOpera.{pdf,png}, compare_Br_z0.png
root -l -b -q 'comparison/compareFieldMaps.C+'          # default data/sphenix3dmapxyz.root; -> plots/01_…–11_…
root -l -b -q 'comparison/checkTilt.C+'                 # -> tilt_A…F_*.pdf (m=1 Bz maps/profiles)
root -l -b -q 'comparison/compareTiltedOpera.C+("data","data/sphenix3dmapxyz.root","plots")'  # -> tiltedOpera_*.{pdf,png}
root -l -b -q 'comparison/findCenterOpera.C+("data/sphenix3dmapxyz.root")'  # console only

# Drop-in Cartesian map for sPHENIX reconstruction
root -l -b -q 'export/makeDeliveredMaps.C+'             # -> output/sphenix_solenoid_{opera_matched_to_mapping,measured_smoothed}_2022-12-02.root
```

Every plot in `plots/` is reproduced by re-running the macro listed above; the
three tilt estimators are intentionally kept as separate macros (see the
[tilt cross-checks](#cross-checks-three-independent-tilt-estimators)):
`checkAlignment.C` (global m=0), `estimateTilt.C` (near-axis + axis-line fit),
and `comparison/checkTilt.C` (m=1 Bz).

The `sPHENIXFieldMap` (r, z) grid is r ∈ [0, 900] mm (25 mm step, 37 nodes),
z ∈ [−2700, 2100] mm (20 mm step, 241 nodes); bilinear interpolation, Bφ = 0,
and Br derived from ∇·B = 0.

## Results

### Measured field map (2022-12-02)

| Quantity | Value |
|----------|-------|
| On-axis field, Bz(0,0,0)     | **1.397 T** |
| Peak on-axis Bz              | 1.3975 T |
| Magnetic centre (Br zero crossing) | **z ≈ +26 mm** |
| ∇·B residual (RMS / \|max\|) | 2.7 × 10⁻⁶ / 2.0 × 10⁻⁵ T/mm |
| Solenoid tilt \|θ\|          | ≈ 4.0–4.7 mrad toward −x (azimuth ≈ 176°) **apparent at full field** — the [half-field cross-check](#half-field-cross-check-a-second-excitation) splits this into ~2.2 mrad field-proportional + a ~1.4 mT current-independent offset (survey ≈ 0.17 mrad); see [tilt cross-checks](#cross-checks-three-independent-tilt-estimators) and [Global synthesis](#global-synthesis) |

### Comparison with OPERA

The official OPERA map (`sphenix3dmapxyz.root`, the `shiftby2p85cm` 28.5 mm
coil-offset default) is `81×81×101` over x,y ∈ [−80, 80] cm, z ∈ [−100, 100] cm.
The comparison window is the tracking volume (r ≤ 80 cm, |z| ≤ 100 cm).

| Quantity | Measured (2022-12-02) | OPERA |
|----------|-----------------------|-------|
| Magnetic centre (z) | ~+26 mm | ~+28.5 mm |
| Bz(0,0,0)           | 1.397 T | 1.385 T |
| Peak Bz             | 1.3975 T | 1.3848 T |

Both maps agree with each other and with the **28.5 mm design coil offset**
(`shiftby2p85cm`); see [Magnetic centre](#magnetic-centre) below for the
per-estimator values.

- **On-axis ΔBz(z=0) = −12.7 mT** (OPERA − measured): the measured central field
  is ~0.9 % higher than OPERA, a nearly uniform scale offset.
- **Max |ΔBz| = 13.9 mT**, **max |ΔBr| = 3.1 mT** over the tracking volume — no
  z-offset or shape disagreement.
- **Azimuthal (m=1) structure is small:** Bz m=1 amplitude ≤ 0.88 mT, with a
  consistent phase (Bφ m=1 phase 5.9 ± 0.8°), i.e. a single rigid-body
  tilt/offset direction rather than a winding asymmetry.
- **Maxwell residuals:** measured |∇·B| RMS = 0.003 mT/cm (enforced by
  construction) vs OPERA 45.9 mT/cm — the latter is the known artifact of storing
  Bx,By,Bz as three independent trilinear grids, not a defect of either map.

Residual RMS of `|measured − s·OPERA|` over the tracking volume (~102 k raw
measured points), from `comparison/compareTiltedOpera.C`:

| corrections applied | scale s | residual RMS \|B\| | what remains |
|---------------------|:-------:|:-----------------:|--------------|
| none (raw difference)       | 1.000  | 16.9 mT | scale + rotation |
| normalization only          | 1.0091 | 11.4 mT | rotation (~6 mT) |
| rotation only               | 1.000  | 15.8 mT | scale (~13 mT) |
| **normalization + rotation**| 1.0091 | **9.8 mT** | transverse noise |

At the final stage the residual splits into **Bz = 0.84 mT** (the meaningful
number — better than one part in a thousand) and **transverse = 9.7 mT**, which is
point-to-point measurement noise: the *coherent* transverse field is the 4.4 mrad
tilt, already removed in azimuthal average (the m=1 amplitude is ≤ 0.88 mT). The
best-fit scale 1.0091 confirms the +0.9 % seen in Bz(0,0,0).

### Magnetic centre

On-axis Bz is flat to < 0.1 mT over several cm around the peak, so a plain
`argmax` on the field grid cannot localise the centre (on the 2 cm OPERA grid the
+20 and +40 mm nodes are equal to 1×10⁻⁶ T, so argmax reports +40 mm). The centre
must be found from the shape of the curve. `analysis/findCenter.C` (measured map)
and `comparison/findCenterOpera.C` (OPERA) report the same three robust estimators
plus the naive argmax for reference:

| Estimator | Measured (2022-12-02) | OPERA |
|-----------|-----------------------|-------|
| Bz argmax (grid-limited) | +20 mm ⚠️ | +40 mm ⚠️ |
| **Bz symmetry**          | **+27.0 mm** | **+31.0 mm** |
| **Parabola vertex**      | **+27.0 mm** | **+30.7 mm** |
| **Br zero crossing** (r = 100 mm) | **+25.8 mm** | **+32.2 mm** |

The three robust estimators agree within a few mm for each map, and both maps land
on the **28.5 mm design coil offset** — confirming the apparent measured-vs-OPERA
centre difference seen with `argmax` (+20 vs +40 mm) was an artifact, not physics.

```bash
root -l -b -q 'analysis/findCenter.C+("data")'
root -l -b -q 'comparison/findCenterOpera.C+("data/sphenix3dmapxyz.root")'
```

### Solenoid tilt: measured vs OPERA

A rigid tilt of the solenoid axis relative to the sPHENIX z axis shows up as a
non-zero **φ-averaged transverse field**: for a perfectly axial field the radial
component cancels in the azimuthal mean, so a residual ⟨B⊥⟩ ≈ |B|·θ measures the
tilt. `analysis/checkAlignment.C` computes this for the measurement; the same
calculation applied to the OPERA `TNtuple` (mean of `bx,by,bz` over the symmetric
Cartesian grid) gives the OPERA tilt.

| Map | \|tilt\| θ | direction | net ⟨B⊥⟩ |
|-----|-----------|-----------|----------|
| **Measured (2022-12-02)** | **4.13 mrad** (θx −4.12, θy +0.29) | azimuth ≈ 176° (toward −x) | ⟨Bx⟩ ≈ −5.8 mT |
| **OPERA** (official `sphenix3dmapxyz.root`) | **≈ 0.02–0.05 mrad** | (chimney only) | ⟨Bx⟩,⟨By⟩ < 0.07 mT |

**The measured solenoid carries a real ~4 mrad axis tilt; OPERA does not.** OPERA
is an idealized calculation, axisymmetric apart from the chimney, so its residual
transverse field is two orders of magnitude smaller (and is just the chimney
breaking azimuthal symmetry). The measured tilt is a physical misalignment of the
magnet axis relative to the surveyor frame and **should be cross-checked against
the survey before being quoted** as a hardware number.

**Physical picture.** With the sPHENIX site convention (+z = north, +x = west,
y = up, so −x = east), θy ≈ 0 means the magnet is **level** (no pitch/roll); the
tilt is a pure horizontal yaw of ~4 mrad about the vertical axis, swinging the
**north end toward the east** (south end toward the west). The field follows the
axis: mostly north, canting slightly east (⟨Bx⟩ ≈ −6 mT eastward on the ~1.4 T
northward field). The ~4 mrad (≈ 0.23°) corresponds to ~6 mm of axis displacement
over a ~1.5 m half-length.

The same tilt appears as a small m=1 modulation of Bz: `comparison/checkTilt.C`
finds a mean m=1 Bz phase of ~7° in the central tracking region with amplitude
< 1 mT, and `compareFieldMaps.C` finds the OPERA−measured difference is dominated
by a single (m=1) direction (Bφ m=1 phase 5.9 ± 0.8°), i.e. a rigid-body
tilt/offset rather than a winding asymmetry.

> **Note — correction to the earlier comparison.** The retired
> `cern-opera-comparison` repo stated that *OPERA* contained a ~4.1 mrad dipole
> "toward −y" that was absent from the measurement. Measured against the **official**
> OPERA map (`sphenixoperamaps/sphenix3dmapxyz.root`) with the φ-averaged method
> above, OPERA shows **no** such tilt — it is the **measurement** that carries the
> ~4 mrad tilt. The old statement most likely referred either to a different,
> derived OPERA file (that analysis was hardwired to a 111³ `…trackingmapxyz.root`
> with an extra `hz` branch, not the official map) or to the m=1 phase of the
> *difference* rather than a tilt of OPERA itself.

#### Cross-checks: three independent tilt estimators

The φ-averaged method above (`checkAlignment.C`) volume-weights toward large
radius and gives a single global number. `analysis/estimateTilt.C` adds two
further estimators that use **physically distinct signatures** of the same rigid
tilt, so they fail in different ways — if all three agree, the tilt is real and
not an artefact of one method. All use only *measured* points (no r = 0
extrapolation through the field-map class).

**Method 1 — near-axis m=0 transverse field.** The φ-averaged transverse field
cancels the axisymmetric radial term and leaves ⟨B⊥⟩ ≈ θ·Bz at every radius, so
θ(r) should be flat for a rigid tilt. Per single ring it is **not** flat — it
scatters from 1.7 to 12 mrad because the ~5 mT transverse signal is comparable to
the ring-to-ring measurement noise (the innermost r = 50 mm ring alone gives a
spurious 6.6 mrad). Averaging the inner rings (r ≤ 300 mm, |z| < 150 mm) gives
**4.65 ± 0.46 mrad toward φ ≈ −171°**, where the error is the standard error of
the mean across the 11 independent rings — i.e. the ring-to-ring scatter is the
dominant *statistical* uncertainty (the larger systematic is discussed below).
This is the direct, least model-dependent estimate, and it confirms that a
*single* small volume at the centre is too noisy to trust.

**Method 2 — magnetic-axis line + tilt fit.** This separates a tilt from a pure
translation of the axis, which the φ-average alone cannot. Near the axis the
transverse field is the sum of a uniform tilt term and the solenoid focusing
field that points back toward the magnetic axis:

```
Bx(x,y,z) ≈ θx·B(z)  −  ½ B′(z)·( x − x_axis(z) ),   x_axis(z) = x0 + θx·z
```

with B′ = dBz/dz (and analogously for y). A per-z linear fit `Bx = px(z) + g(z)·x`
across the bore measures the focusing gradient `g(z) = −½B′(z)` and the intercept
`px(z)` at each z slab — driven by the fringe, where g is large, using real points
across the measured bore rather than an extrapolation to r = 0. A global
least-squares of the intercepts,

```
px(z) = θx·( B(z) − g(z)·z )  −  x0·g(z),
```

then solves for the tilt slope θx and the axis offset x0 simultaneously (likewise
θy, y0). Errors come from the least-squares covariance (residual variance × the
parameter covariance matrix). The fit gives **|θ| = 4.40 ± 0.08 mrad toward
φ ≈ −178°** (θx = −4.40 ± 0.08, θy = −0.16 ± 0.03 mrad) **plus a small transverse
axis offset (x0, y0) = (+2.6 ± 0.5, +1.2 ± 0.2) mm** at z = 0, i.e. the tilt and
the residual translation are disentangled.

Each estimator carries a statistical (internal-precision) error only:

| Estimator | \|θ\| [mrad] | direction φ | stat. error source (internal) |
|-----------|--------------|-------------|-------------------------------|
| `checkAlignment.C` — global φ-averaged m=0 | 3.99 ± 0.09 (4.13 point-weighted) | 176.2 ± 1.6° | spread of strong-field z-slices (\|Bz\|>1 T) |
| `estimateTilt.C` M1 — near-axis, mean of 11 rings (r ≤ 300 mm) | 4.65 ± 0.46 | −171 ± 6° | ring-to-ring standard error |
| `estimateTilt.C` M2 — axis-line fit | 4.40 ± 0.08 | −178.0 ± 0.3° | least-squares fit covariance |

All three land at **4.0–4.7 mrad in the −x direction** (θx ≈ −4.4 mrad, θy ≈ 0).

**The ± values above are statistical only and should not be read as the
uncertainty on a hardware tilt.** All three are *dispersion*-based (scatter across
z-slices, across rings, or about the fit model), so they measure only how
reproducibly *this* map pins down the field-to-frame rotation. The dominant
systematic — a probe/mapper/registration **yaw of the measurement frame** — is
**common-mode**: it adds the same offset to every slice, ring and z-slab, so it
contributes *exactly zero* to any of these dispersions. The error bars are
therefore structurally **blind** to the effect most likely to fake the signal;
they would read ~0.1 mrad even if a 4 mrad frame yaw produced the entire result
(field data cannot separate a magnet yaw from a measurement-frame yaw — see
[Systematics](#systematics-is-the-4-mrad-tilt-real)).

The realistic uncertainty hierarchy is:

| source | size | what it captures |
|--------|------|------------------|
| per-method statistical (above) | 0.1–0.5 mrad | random precision on this map |
| method-to-method spread | ~0.5 mrad | estimator dependence (3.99 / 4.40 / 4.65) |
| frame-yaw degeneracy | up to ~100% | cannot be bounded from the map alone |

So the honest statement is: **measured field-to-frame angle ≈ 4.4 mrad (stat. ~0.1
mrad), with a systematic uncertainty of order the value itself — consistent with
the magnet being aligned to the survey's ~0.2 mrad.** The map establishes that a
~4 mrad rotation exists between field and nominal frame; it cannot by itself say
whether that rotation is the magnet, the probe, or the registration. The survey is
the decisive external check. (Of the three statistical errors, M1's ±0.46 is the
most trustworthy — ring-to-ring is genuinely semi-independent; M2's ±0.08 is
optimistic even statistically, since the p(z) residuals carry unmodelled
structure.) See `plots/tilt_estimates.{pdf,png}`.

### Systematics: is the ~4 mrad tilt real?

The estimators agree that a ~4 mrad rotation exists between the measured field and
the nominal (surveyor) frame. What they **cannot** establish is *what* is rotated.

**The fundamental degeneracy.** A rigid rotation of the *measuring system*
relative to the survey frame produces field data identical to rotating the
*magnet* by the same angle — `B_measured = R · B_true` looks the same whether `R`
came from the magnet yawing one way or the probe/mapper yawing the other. From the
field map alone, "magnet yawed east" and "probe yawed west" are indistinguishable.
The signal is also small in fractional terms: ⟨B⊥⟩ ≈ 6 mT on a 1.4 T field is
**0.4 %**, i.e. 4 mrad × 1.4 T — so a few-mrad error *anywhere* in the chain
reproduces it exactly.

**Candidate effects.**

| effect | ∝ Bz? | fakes a yaw? | notes |
|--------|:---:|:---:|-------|
| Probe-triad mounting rotation | yes | yes | mechanical misalignment leaks Bz into B⊥ |
| Mapper/gantry yaw vs. fiducials | yes | yes | whole field appears rotated |
| Fiducial / registration error | yes | yes | old map already had a 265 mm *z* survey error |
| Rotation-stage tilt / readout offset | yes | yes | if the probe is spun in φ to sample azimuth |
| Probe transverse–axial cross-talk | yes | yes | planar Hall effect; ~4 mrad ≈ whole signal |
| Mapping-arm deflection (magnetic force) | yes | yes | field-dependent; absent in a magnet-off survey |
| External / ambient uniform field | **no** | no | disfavoured: ⟨B⊥⟩ tracks Bz |
| Coil winding asymmetry | yes | no *(real field)* | built-in transverse dipole, not a yaw; different z-dependence (`checkTilt.C`) |

**What points toward a measurement systematic rather than a hardware yaw:**

- **The survey says ~0.17 mrad** mechanically (pending the new one) — two orders of
  magnitude smaller.
- **It is purely horizontal** (θy ≈ 0). Gravity-driven effects (cryostat/coil sag)
  would appear as a *vertical* pitch, not a horizontal yaw; a horizontal-only
  rotation is more naturally an installation/registration error about the vertical
  axis or a mapper yaw.
- **⟨B⊥⟩ ∝ Bz** through the plateau (θ roughly constant in z) — the signature of a
  rotation/cross-talk, which is what rules out an additive external field.

**Diagnostics that could actually disentangle it:**

- **The pending cryostat survey** — the only external handle that breaks the
  degeneracy. Decisive.
- **Current dependence** — if maps at different excitations exist: does ⟨B⊥⟩ scale
  linearly through zero with Bz (→ rotation / cross-talk / deflection) or carry a
  current-independent offset (→ external field)?
- **Reproducibility across a deliberate remount** — re-survey, re-mount, re-map; if
  the yaw changes, it lives in the measurement chain.
- **Probe calibration residual** — the transverse reading in a known pure-axial
  field gives the cross-talk floor directly.

> **Note.** Comparing the fine and rough scans is **not** an independent
> cross-check: the rough map was taken first in the *same campaign and setup* as
> insurance, so it shares the same probe alignment, mapper registration and any
> field-dependent deflection. Both inherit an identical frame yaw.

**Bottom line.** Treat the ~4 mrad as a rotation between field and frame whose most
likely home is the measurement chain. The magnet-yaws-east interpretation can be
neither confirmed nor excluded from the map alone; the survey is required.

### Half-field cross-check: a second excitation

The full-field fine and rough scans cannot test the tilt (they share the same
setup and frame yaw — see the Note above). The one lever they lack is a
**different excitation**, and the 2022-12-02 campaign recorded one: a **half-field**
rough map (`pointCloudRoughHalfField.csv`) on the *identical* 10 cm / 10° grid as
the full-field rough map — same probe, mapper and registration, only the current
changed (plateau ratio verified **0.500**). `analysis/checkHalfField.C` runs the
comparison; the full-field analysis is untouched.

**Why a second excitation is decisive in principle.** A genuine rigid tilt is a
rotation `B_meas = R·B_true`, so the transverse leakage is *exactly* proportional
to Bz and the angle θ = ⟨B⊥⟩/Bz is the **same at any current**. Anything whose
transverse field does **not** scale with Bz is, by definition, not a rigid tilt.

Three observations, all consistent:

1. **The apparent tilt is not invariant.** The global φ-averaged angle **grows**
   from **3.09 mrad at full field to 4.14 mrad at half field**, and ⟨B⊥⟩ falls
   only to **0.67** of its full-field value — not the **0.50** a pure rotation
   requires. (The axis-line fit agrees: θ_x −3.48 → −4.62 mrad.)
2. **The magnetic centre does not move.** Between the two excitations the
   transverse axis (x₀, y₀) shifts by (+0.5, −0.4) mm — within the ~1 mm fit
   error — and the axial centre by **0.2 mm**. An additive transverse offset
   biases the *fitted tilt* (it projects onto the ∝Bz term) but not the
   focusing-derived centre, so a **growing angle with a stationary centre** is its
   fingerprint.
3. **A through-zero fit separates the two pieces.** Pooling the flat-field slices
   (|z| < 700 mm, where the focusing term vanishes) from both excitations gives a
   clean line that **misses the origin**:

   - ⟨Bx⟩ = (**−2.14 ± 0.05 mrad**)·Bz + (**−1.39 ± 0.05 mT**)
   - ⟨By⟩ = (**−0.65 ± 0.07 mrad**)·Bz + (**+0.28 ± 0.07 mT**)

   The slope is a **field-proportional yaw of 2.24 mrad** (φ ≈ −163°); the
   intercept is a **current-independent transverse offset of 1.42 mT** (φ ≈ +169°),
   **nonzero at 26σ**. The 0.09 mT residual rms shows the linear model is excellent
   (no curvature → no significant iron-saturation nonlinearity).

![Half-field through-zero fit](plots/halffield_BperpVsBz.png)

**What it establishes.** The full-field "≈4 mrad tilt" is the sum of two physically
distinct things: a **field-proportional yaw of ~2.2 mrad** and a
**current-independent ~1.4 mT transverse offset** that the full-field map alone
reported as extra tilt. This **revises** the earlier "external/ambient field
disfavoured" entry in the candidate table — that was inferred from the
z-flatness of θ *within one excitation*; the half-field map measures the
current-independent term directly and finds it nonzero. Such a term is most
naturally a **transverse Hall-probe pedestal** (a zero offset); it cannot be
distinguished from a real remanent/ambient ~1.4 mT field by the map alone — a
**magnet-off / zero-field probe reading**, which this campaign did **not** take,
is what would settle it. **What it does not do:** the 2.24 mrad slope still mixes a
real magnet/frame rotation with any ∝Bz probe cross-talk, so even that residual is
an upper bound on a physical yaw.

### Global synthesis

Combining the full-field estimators, the half-field decomposition and the survey,
the transverse-field "tilt" resolves into three layers of very different status:

| component | size | nature | scales as |
|-----------|------|--------|:---------:|
| current-independent offset | ~1.4 mT (~1.0 mrad-equiv at full field) | probe pedestal, or remanent/ambient field | constant |
| field-proportional yaw | 2.24 ± 0.05 mrad | rigid rotation **and/or** ∝Bz probe cross-talk | ∝ Bz |
| real magnet yaw (survey) | ≈ 0.17 mrad | mechanical alignment | geometric |

Read top-to-bottom these are successive subtractions: raw apparent ~4 mrad →
remove the offset → **2.24 mrad** field-proportional → remove cross-talk (size
unknown) → the survey's **~0.17 mrad**. So the honest statement is:

> **The measured field shows a ~4 mrad apparent transverse tilt at full
> excitation, of which only ~2.2 mrad scales as a rotation and the rest is a
> current-independent ~1.4 mT offset (most likely a probe pedestal). Whether even
> the 2.2 mrad is a real magnet yaw or further measurement systematic cannot be
> decided from the maps; the cryostat survey (~0.17 mrad) is the external check.**

**On the errors.** The per-method dispersion errors quoted above (±0.08–0.46 mrad)
measure only how reproducibly each estimator pins the angle *on one map*; they are
structurally blind to the common-mode systematics that dominate. The half-field
result makes that concrete — the *same magnet* read 3.1 vs 4.1 mrad just by halving
the current. The realistic uncertainty on any "magnet tilt" is therefore **of order
the value itself**, not the printed ±.

**Consequence for the delivered map.** Because the maps cannot collapse this to one
number, the delivery is a **menu over the yaw**, all sharing the +0.9 % amplitude
(built by `export/makeDeliveredMaps.C+(yawMrad)`; file-level details and md5s in
[`export/MANIFEST_fieldmaps.txt`](export/MANIFEST_fieldmaps.txt)):

| yaw | rationale |
|----:|-----------|
| **0.0 mrad** | no-tilt / scale-only — consistent with the survey and with the apparent tilt being a measurement systematic |
| **2.2 mrad** | field-proportional yaw, pedestal removed — best single estimate **if** a real rotation is assumed (an upper bound) |
| **4.4 mrad** | raw apparent yaw — **disfavoured**; ~1 mrad of it is the additive offset, and the remainder still exceeds the survey |

The **0.0 and 2.2 mrad maps bracket the defensible range**; 4.4 mrad is retained
only for continuity with the earlier single-yaw delivery. A survey would collapse
the menu to one.

### Maps for reconstruction

`export/makeDeliveredMaps.C+(yawMrad)` writes the drop-in maps described in
[Maps delivered to tracking](#maps-delivered-to-tracking) — all `TNtuple`
(`x:y:z:bx:by:bz:hz`, cm/T) on the production 111³ grid (±110 cm, 2 cm steps)
used by `PHField3DCartesian`, carrying the measured +0.9 % amplitude and a chosen
horizontal yaw. Run it once per yaw to produce the menu (0.0, 2.2, 4.4 mrad); the
default reproduces the original 4.4 mrad maps exactly. The deployment note that
ships alongside the maps is
[`export/MANIFEST_fieldmaps.txt`](export/MANIFEST_fieldmaps.txt) (the six md5s, the
yaw menu, and the flag that the currently-deployed CDB map is z-shifted — centre at
z ≈ −240 mm — and must be replaced).

## Summary

> ### The measured 2022-12-02 field validates the OPERA map used in reconstruction.
>
> After a single uniform **0.9 % normalization**, the measured **Bz** reproduces
> the OPERA calculation to **0.84 mT RMS** point-by-point across the tracking
> volume — **better than one part in a thousand**, with no z-offset and no shape
> disagreement (see [Comparison with OPERA](#comparison-with-opera)). A few-mrad
> transverse yaw changes Bz only at the µT level, so this validation is independent
> of the tilt question.
>
> The physical differences are therefore **one calibration number** and **one
> alignment question**:
>
> - **~0.9 % scale** — the measured field is ~0.9 % above OPERA; a probe-calibration
>   question for the CERN mapping group.
> - **a transverse tilt** — a ~4 mrad apparent yaw at full field, which the
>   [half-field cross-check](#half-field-cross-check-a-second-excitation) resolves
>   into a **~2.2 mrad field-proportional** part plus a **current-independent
>   ~1.4 mT offset** (most likely a probe pedestal); the survey says ~0.17 mrad. The
>   delivered maps therefore span a **yaw menu (0 / 2.2 / 4.4 mrad)** rather than a
>   single angle (see [Global synthesis](#global-synthesis)).

## Maps delivered to tracking

All delivered maps are on the production 111³ / ±110 cm `PHField3DCartesian` grid
(`x:y:z:bx:by:bz:hz`, cm/T) and carry the **+0.9 % amplitude scale** (θ_y ≈ 0
throughout, so the magnet is treated as level). They come in two **provenance
variants** × the three-point **yaw menu**, built by
`export/makeDeliveredMaps.C+(yawMrad)`. File-level details and md5s are in
[`export/MANIFEST_fieldmaps.txt`](export/MANIFEST_fieldmaps.txt).

**Provenance variant** — pick by coverage need:

- **`opera_matched_to_mapping`** — the OPERA calculation **rescaled ×1.0091** and
  yawed to match the mapping; corrected in the spirit of the measurement, **not** a
  measured map. **Defined everywhere** in the cube — use it for anything sampling
  r > 90 cm (field integrals, swimming to calorimeter, TPC outer edge).
- **`measured_smoothed`** — the **measured** field: the φ-averaged, ∇·B-enforced
  (r,z) reconstruction (real amplitude and z-profile, ~10 mT per-point transverse
  noise removed) with the yaw reinserted as a rigid rotation. **Zero beyond
  r = 90 cm**, where there is no measurement.

The two agree to **~0.84 mT** in the tracking volume at a given yaw (that
near-identity is the [validation result](#comparison-with-opera)); they differ only
in provenance and in how the unmeasured r > 90 cm corners are treated (OPERA-filled
vs. zero).

**Yaw menu** — each value is delivered as *both* provenance variants, named
`sphenix_solenoid_{opera_matched_to_mapping,measured_smoothed}_2022-12-02_yawX.Xmrad.root`
(see [Global synthesis](#global-synthesis) for why a menu):

| yaw | file suffix | rationale |
|----:|-------------|-----------|
| 0.0 mrad | `_yaw0.0mrad` | no-tilt / scale-only — the bare amplitude match; consistent with the survey |
| 2.2 mrad | `_yaw2.2mrad` | field-proportional yaw, pedestal removed — best single estimate **if** a real rotation is assumed |
| 4.4 mrad | `_yaw4.4mrad` | raw apparent yaw — **disfavoured** (~1 mrad is the additive offset); kept for continuity |

Which to register for production is a tracking-group decision. The **0.0 and
2.2 mrad maps bracket the defensible range**; 4.4 mrad is retained only for
continuity with the earlier single-yaw delivery.

Notes:

- All six carry the +0.9 % scale baked in, so run them with **`magfield_rescale =
  1.0`**. (The yaw-0.0 `opera_matched` map *is* the old "scale-only" recipe — bare
  OPERA ×1.0091 — now shipped as a file.)
- The **measured map deliberately omits nothing it measured but keeps no noise**:
  it is the measurement smoothed, *with* its yaw. It is **not** the raw point cloud
  (whose ~10 mT transverse scatter is measurement noise, not field).
- The yaw is applied as a rigid rotation; rebuilding any case is one command. A
  magnet-off/zero-field probe reading or the cryostat survey would collapse the
  menu to one map (see [Global synthesis](#global-synthesis)).

This supersedes the earlier single-yaw delivery and the axisymmetric
`sphenix_measured_fieldmap_cartesian.root` (no tilt, cut off at r = 90 cm); both are
retired.

## Plots

### Measured field overview
![Field overview](plots/fieldMap_overview.png)

### Measured vs OPERA (overlays)
![Measured vs OPERA](plots/compare_newVsOpera.png)

### Solenoid alignment (global m=0 transverse field vs z)
![Alignment](plots/fieldMap_alignment.png)

### Solenoid tilt estimators (near-axis + axis-line fit)
![Tilt estimates](plots/tilt_estimates.png)

### Yawing OPERA by the measured tilt
![Tilt vs z, OPERA yawed to match measured](plots/tiltedOpera_transverse_vs_z.png)
![Field the yaw adds to OPERA](plots/tiltedOpera_yaw_correction_maps.png)

*`comparison/compareTiltedOpera.C` applies the measured ≈ 4.4 mrad yaw
(`R_y(θ)`, θx = −4.40 mrad, θy ≈ 0) to the OPERA map and asks what the
measured − calculated difference would look like. Top: the φ-averaged transverse
field ⟨Bx⟩(z). The measured raw data (black) sit at ≈ −6 mT; OPERA untilted (blue)
is ≈ 0; OPERA yawed by 4.4 mrad (red) reproduces the measured value across the
plateau — the transverse difference is the tilt. Bottom: the field the yaw adds to
OPERA — a uniform ΔBx ≈ −6 mT (left, = the difference it removes) and only a
sub-few-mT, antisymmetric-in-x ΔBz in the fringe corners (right). The ~0.9 %
(~12.7 mT) Bz scale offset is rotation-invariant and is left untouched, so after the
yaw the residual is just that scale offset plus sub-mT fringe structure.*

### Half-field cross-check (tilt vs excitation)
![Half-field tilt and transverse field vs z](plots/halffield_tilt.png)

*`analysis/checkHalfField.C` overlays the full- (blue) and half-field (red)
campaigns. Left: θ_x(z); right: ⟨|B⊥|⟩(z). A rigid tilt would give the same θ at
both excitations and ⟨|B⊥|⟩ scaling by exactly ½; instead θ grows at half field and
⟨|B⊥|⟩ falls by less than ½ — the signature of a current-independent offset on top
of the field-proportional yaw. The through-zero fit that separates the two
(plots/halffield_BperpVsBz.png) is shown in the
[half-field cross-check](#half-field-cross-check-a-second-excitation) section.*

The `comparison/compareFieldMaps.C` macro additionally writes the numbered series
`plots/01_…`–`plots/11_…` (on-axis Bz, 2-D maps, ΔBz, Br comparison, m=1
amplitude/phase, Bφ, radial profiles, φ-dependence, ΔBz(φ=0 vs 180), Maxwell
residuals), and `comparison/checkTilt.C` writes the `plots/tilt_A…F_*.pdf` m=1
Bz maps and profiles.

## Dependencies

- [ROOT](https://root.cern) (developed with 6.40)
- The two large input files above (not distributed in git).
