"""See swell in SWOT: band-passed 250 m Unsmoothed SSHA and its 2D wavenumber spectrum."""
import earthaccess
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

PASS = "059_164"           # cycle_pass, 18 Nov 2026 ~05 UTC, ~137W (from swell_finder.py + tile scan)
DATES = ("2026-11-18", "2026-11-19")
LAT_RANGE = (35.5, 39.5)   # where the swell peak shows up in the spectra
TILE_LAT = 37.8            # centre latitude of the zoom tile
TILE_KM = 40               # tile size for the spectrum
ZOOM_KM = 15               # zoom shown in the middle panel
BAND = (1.2, 2.0)          # search for the peak at |k| in cycles/km, i.e. wavelengths 500-830 m
                           # (KaRIn noise sits mostly at |k| < 1.2; ~2.0 is the Nyquist limit of the pixel grid)

earthaccess.login()
results = earthaccess.search_data(short_name="SWOT_L2_LR_SSH_D", granule_name=f"*Unsmoothed_{PASS}*", temporal=DATES)
f = earthaccess.download(results, "data/")[0]

ds = xr.open_dataset(f, group="right")  # Unsmoothed files hold each swath in its own group ("left", "right")
lat = ds.latitude.median("num_pixels").values
x = ds.cross_track_distance.median("num_lines").values / 1000  # km from nadir (edge pixels have no geolocation)
cols = np.where(np.isfinite(x) & (x > 8))[0]                    # skip pixels nearest nadir
ssha = (ds.ssha_karin_2 + ds.height_cor_xover).where(ds.ssha_karin_2_qual < 2**17)

# segment for context
seg = (lat > LAT_RANGE[0]) & (lat < LAT_RANGE[1])
s_seg = ssha[seg][:, cols].values
y_seg = lat[seg]

# pixel spacing, measured rather than assumed: cross-track is 250 m, but along-track is only ~237 m
i0 = int(np.argmin(abs(lat - TILE_LAT)))
DX = float(np.median(np.diff(x[cols])))                       # km
la, lo = np.radians(ds.latitude[i0 - 50:i0 + 50, cols[0]]), np.radians(ds.longitude[i0 - 50:i0 + 50, cols[0]])
DY = float(np.median(6371 * np.hypot(np.diff(la), np.diff(lo) * np.cos(la[1:]))))  # km (local flat-earth)
print(f"pixel spacing: cross-track {DX * 1000:.1f} m, along-track {DY * 1000:.1f} m")

# tile of TILE_KM x TILE_KM (so more rows than columns); remove scales > ~10 km (eddies etc.); NaNs -> 0
nx, ny = int(TILE_KM / DX), int(TILE_KM / DY)
r0 = i0 - ny // 2
tile = ssha[r0:r0 + ny, cols[:nx]].values
tile = np.nan_to_num(tile - np.nanmean(tile))
tile = tile - gaussian_filter(tile, 8)

# 2D FFT: separate wavenumber axes (cycles/km) for the two directions; Hann window reduces edge leakage
kx, ky = np.fft.fftfreq(nx, d=DX), np.fft.fftfreq(ny, d=DY)
KX, KY = np.meshgrid(kx, ky)  # shape (ny, nx), like the tile
kmag = np.hypot(KX, KY)
F = np.fft.fft2(tile * np.outer(np.hanning(ny), np.hanning(nx)))
spec = gaussian_filter(np.fft.fftshift(np.abs(F) ** 2), 1.5)

# spectral peak inside the band -> wavelength -> period via deep-water dispersion L = g T^2 / 2pi = 1.56 T^2
kxs, kys = np.fft.fftshift(kx), np.fft.fftshift(ky)
kk = np.fft.fftshift(kmag)
iy, ix = np.unravel_index(np.argmax(np.where((kk > BAND[0]) & (kk < BAND[1]), spec, 0)), spec.shape)
L = 1000 / kk[iy, ix]
# keep only wavenumbers within 0.3 cycles/km of the peak (and its mirror): the coherent swell component.
# (the evidence is the spectral peak itself; this filtered image just shows what that peak looks like)
near = (np.hypot(KX - kxs[ix], KY - kys[iy]) < 0.3) | (np.hypot(KX + kxs[ix], KY + kys[iy]) < 0.3)
# a ~650 m wave is only ~2.7 pixels long, so draw it on a 4x finer grid: zero-pad the (band-limited) spectrum,
# keeping k=0 at the centre index that ifftshift expects for the 4x larger array
pad = [(2 * m - m // 2, 4 * m - m - (2 * m - m // 2)) for m in (ny, nx)]
G = np.pad(np.fft.fftshift(np.fft.fft2(tile) * near), pad)
swell = np.real(np.fft.ifft2(np.fft.ifftshift(G))) * 16
print(f"peak: k = ({kxs[ix]:.2f}, {kys[iy]:.2f}) cycles/km -> wavelength {L:.0f} m -> period {np.sqrt(L / 1.56):.1f} s")

fig, axes = plt.subplots(1, 3, figsize=(17, 7), gridspec_kw={"width_ratios": [0.6, 1, 1]})
im = axes[0].pcolormesh(x[cols], y_seg, s_seg, cmap="RdBu_r", vmin=-0.3, vmax=0.3, shading="auto")
axes[0].axhspan(lat[r0], lat[r0 + ny], fill=False, lw=1.5)
axes[0].set(xlabel="cross-track (km)", ylabel="latitude", title="SSHA (m), right swath")
plt.colorbar(im, ax=axes[0], shrink=0.6)

zx, zy = int(ZOOM_KM / (DX / 4)), int(ZOOM_KM / (DY / 4))  # show only the centre of the (4x upsampled) tile
cx, cy = (4 * nx - zx) // 2, (4 * ny - zy) // 2
im = axes[1].imshow(swell[cy:cy + zy, cx:cx + zx], cmap="RdBu_r", vmin=-0.03, vmax=0.03, origin="lower",
                    extent=[0, ZOOM_KM, 0, ZOOM_KM])
axes[1].set(xlabel="cross-track (km)", ylabel="along-track (km)",
            title=f"SSHA filtered around the spectral peak, {ZOOM_KM} km zoom at {TILE_LAT}N")
plt.colorbar(im, ax=axes[1], shrink=0.6, label="m")

axes[2].pcolormesh(kxs, kys, spec, cmap="viridis", shading="auto")
for wl in [0.5, 0.6, 1.0]:
    axes[2].add_patch(plt.Circle((0, 0), 1 / wl, fill=False, color="w", lw=0.6))
    axes[2].text(0.7 / wl, 0.7 / wl, f"{wl*1000:.0f} m", color="w", fontsize=8)
axes[2].plot(kxs[ix], kys[iy], "r+", ms=12)
axes[2].set(xlim=(-2, 2), ylim=(-2, 2), aspect="equal", xlabel="k cross-track (cycles/km)",
            ylabel="k along-track (cycles/km)", title=f"2D spectrum: peak λ ≈ {L:.0f} m, T ≈ {np.sqrt(L/1.56):.1f} s")

fig.suptitle(f"SWOT KaRIn 250 m Unsmoothed SSHA, pass {PASS}, {results[0]['umm']['GranuleUR'][36:51]} UTC "
             f"(~{float(ds.longitude[r0, cols[0]]) - 360:.0f}°E)")
plt.savefig(f"figures/swell_{PASS}.png", dpi=150, bbox_inches="tight")
plt.show()
