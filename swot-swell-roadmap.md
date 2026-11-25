# Seeing Swell from Space: a SWOT Roadmap

**Goal:** reproduce the swell "fingerprint" figures from Ardhuin et al. (2024), *GRL*,
"Phase-resolved swells across ocean basins in SWOT altimetry data" (doi:10.1029/2024GL109658),
then build open tooling around it. Each phase ends with a blog post and code on GitHub.

**Key fact:** SWOT's 250 m Unsmoothed SSHA product resolves swell with wavelengths of ~500 m to 1.3 km
(periods ≳ 18 s) at amplitudes down to ~3 cm. Shorter swell is below the resolution floor.

---

## Phase 0 — Setup (one weekend)
- Earthdata login → PO.DAAC access (`earthaccess` in Python); optional AVISO account for L3 products.
- Environment: xarray, dask, numpy, scipy, matplotlib, cartopy.
- Read the L2_LR_SSH product handbook; note the quality flags and the `Unsmoothed` vs `Basic/Expert` distinction.
- **Deliverable:** one SWOT granule opened and plotted.

## Phase 1 — See a swell (2–3 weekends)
- Pick a documented event: ex-Hurricane Kirk swell, NE Atlantic, ~9 Oct 2024, or ex-Hurricane Erin, Bay of Biscay, ~26 Aug 2025.
- Download Unsmoothed L2 (or AVISO L3) SSHA passes over the North Atlantic for those dates.
- Band-pass SSHA between 250 m and 5 km. Swell appears as near-parallel lines; crossing swells criss-cross.
- **Deliverable:** the picture + the ~40 lines that made it. *Blog post 1.*

## Phase 2 — Measure it (1–2 months)
- Tile the swath (~10 km tiles), 2D FFT → directional wavenumber spectrum.
- Peak → wavelength and direction (180° ambiguity resolved from storm position / model).
- Period via deep-water dispersion: λ ≈ 1.56 T² (m, s).
- Compare against WaveWatch III and a buoy the swell reached (Candhis, UK K-buoys, NDBC/CDIP).
- Answer key: AVISO **L3_LR_WIND_WAVE** product (corrected spectrum, swell partition, KaRIn vs WW3 parameters).
- **Deliverable:** `swot-swell` package: granule in → swell spectra out. *Blog post 2:* how close did I get.

## Phase 3 — Reproduce the paper figure (1 month)
- Rebuild the Ardhuin et al. basin-scale figure: multiple passes through one storm's swell field, with wavelength/direction/amplitude per pass.
- Add the instrument transfer-function correction (see Bohé et al. 2025, *IEEE*) so amplitudes match.
- Show dispersion in action: the longest periods outrun the rest; direction fans; amplitude decays.
- **Deliverable:** side-by-side of your figure and the published one. *Blog post 3.*

## Phase 4 — The surfer's question (ongoing)
- Follow one storm's swell with every SWOT pass that crossed it over several days.
- Predict coastal arrival with group velocity c_g ≈ 0.78 T m/s; check against buoy first-arrival.
- Realistic precision at basin scale: ~±3 h, ~±150 km (period ±2.5%, direction ±3°).
- Add Sentinel-1 wave-mode and CFOSAT SWIM (both near-real-time) + WW3 as a live base layer; SWOT backfills ~3 days later.
- **Deliverable:** "I watched this swell cross the Atlantic from orbit." *Blog post 4.*

---

## Moonshots

**M1 — 2D wave-group geometry.** Groups are finite along-crest as well as along-travel; the theory (Longuet-Higgins, 1957/84) has barely been tested. Hilbert-transform the swell field to get the envelope, measure group length/width/run-length statistics across many passes, compare with Gaussian predictions for the observed spectral width.

**M2 — Do sets survive?** Linear theory says a specific group loses identity after ~N² wavelengths (tens of km). Find passes 30–80 km offshore of a buoy that reports raw displacement (CDIP, some Candhis stations). Convert the swath envelope along the propagation direction into a predicted time series; cross-correlate with the buoy 1–2 h later. Linear prediction: no correlation. A peak at the predicted lag = evidence of nonlinear group persistence.

**M3 — Coastal reflection & island diffraction catalogue.** Ardhuin et al. flag reflected swell as visible; nobody has systematically catalogued it. Automated detection of crossing patterns near coasts and archipelagos.

**M4 — Lightweight swell-fingerprint viewer.** Web app rendering per-pass band-passed SSHA over a map, with SWIM/S-1/WW3 overlays. Current tooling is research-grade Python; a browsable viewer would be used.

---

## Resources
- PO.DAAC SWOT L2_LR_SSH (forward-processed `PIC0`, ~3-day latency; `PGC0` reprocessed)
- AVISO L3_LR_SSH and L3_LR_WIND_WAVE + user handbook (THREDDS catalog, free account)
- SWOT Science Team meeting slides: swotst.aviso.altimetry.fr (2025 "Wave Spectra from SWOT")
- LOPS/IFREMER (Ardhuin group) pages and AVISO "Image of the Month" write-ups for Kirk and Erin
- Reference papers: Ardhuin et al. 2024 GRL; Bohé et al. 2025 IEEE; Archer et al. 2025 Nature (for context)

## Budget
Claude Pro ($20/mo, includes Claude Science) · data free · laptop with 16 GB RAM and ~100 GB free; optional 1 TB SSD (~$70) · cloud compute optional (~$20–40 if used, us-west-2 to avoid egress).
