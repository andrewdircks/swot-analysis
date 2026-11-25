"""Scan one SWOT pass for swell: 2D spectra of 40 km tiles at several latitudes, both swaths.

Swell = a compact spectral peak outside KaRIn's noise blob (|k| > ~1.2 cycles/km) that shows up at the SAME
place in neighbouring tiles and on BOTH swaths. Noise only gives random speckles that move from tile to tile.
"""
import earthaccess
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

PASS = "059_164"
DATES = ("2026-11-18", "2026-11-19")
LATS = [33, 34.5, 35.5, 36.3, 37, 37.8, 38.6, 39.5, 41, 43, 45, 47]  # tile centres
TILE_KM = 40
BAND = (1.2, 2.0)  # peak search range, cycles/km (outside the noise blob, inside Nyquist)

earthaccess.login()
results = earthaccess.search_data(short_name="SWOT_L2_LR_SSH_D", granule_name=f"*Unsmoothed_{PASS}*", temporal=DATES)
f = earthaccess.download(results, "data/")[0]

fig, axes = plt.subplots(2, len(LATS), figsize=(2.2 * len(LATS), 5), layout="constrained")
print(" side   lat    lon     std(m)  SNR   k_cross  k_along  wavelength  period")
for row, side in enumerate(["left", "right"]):
    ds = xr.open_dataset(f, group=side)
    lat = ds.latitude.median("num_pixels").values
    x = ds.cross_track_distance.median("num_lines").values / 1000  # km from nadir (left < 0 < right)
    # pixels 8-48 km from nadir, sorted by x so both swaths have the same orientation (no mirrored spectra)
    cols = np.where(np.isfinite(x) & (abs(x) > 8))[0]
    cols = cols[np.argsort(abs(x[cols]))]
    DX = float(np.median(abs(np.diff(x[cols]))))
    nx = int(TILE_KM / DX)
    cols = cols[:nx][np.argsort(x[cols[:nx]])]

    for j, tile_lat in enumerate(LATS):
        ax = axes[row, j]
        i0 = int(np.nanargmin(abs(lat - tile_lat)))
        # along-track spacing measured from lat/lon (~237 m, not 250 m)
        la, lo = np.radians(ds.latitude[i0 - 50:i0 + 50, cols[0]]), np.radians(ds.longitude[i0 - 50:i0 + 50, cols[0]])
        DY = float(np.median(6371 * np.hypot(np.diff(la), np.diff(lo) * np.cos(la[1:]))))
        ny = int(TILE_KM / DY)
        r0 = i0 - ny // 2

        tile = ds.ssha_karin_2[r0:r0 + ny, cols].where(ds.ssha_karin_2_qual[r0:r0 + ny, cols] < 2**17).values
        if np.isnan(tile).mean() > 0.1:  # too much missing data (land, flags)
            ax.set_axis_off()
            continue
        tile = np.nan_to_num(tile - np.nanmean(tile))
        tile = tile - gaussian_filter(tile, 8)  # remove scales > ~10 km

        kx, ky = np.fft.fftshift(np.fft.fftfreq(nx, DX)), np.fft.fftshift(np.fft.fftfreq(ny, DY))
        kk = np.hypot(*np.meshgrid(kx, ky))
        spec = np.fft.fftshift(np.abs(np.fft.fft2(tile * np.outer(np.hanning(ny), np.hanning(nx)))) ** 2)
        spec = gaussian_filter(spec, 1.5)

        # peak in the band; SNR = peak power / median power on the same ring (|k| within 0.1 of the peak's),
        # i.e. how much the peak stands out from everything else at that wavelength
        iy, ix = np.unravel_index(np.argmax(np.where((kk > BAND[0]) & (kk < BAND[1]), spec, 0)), spec.shape)
        ring = abs(kk - kk[iy, ix]) < 0.1
        snr = spec[iy, ix] / np.median(spec[ring])
        L = 1000 / kk[iy, ix]
        print(f"{side:5s} {tile_lat:5.1f}  {float(ds.longitude[i0, cols[0]]) - 360:6.1f}  {np.std(tile):6.3f} "
              f"{snr:5.1f}  {kx[ix]:6.2f}   {ky[iy]:6.2f}    {L:5.0f} m    {np.sqrt(L / 1.56):4.1f} s")

        ax.pcolormesh(kx, ky, spec, cmap="viridis", shading="auto")
        for wl in [0.5, 1.0]:
            ax.add_patch(plt.Circle((0, 0), 1 / wl, fill=False, color="w", lw=0.5))
        ax.set(xlim=(-2, 2), ylim=(-2, 2), aspect="equal", xticks=[-2, 0, 2], yticks=[-2, 0, 2])
        ax.set_title(f"{side} {tile_lat}N\nSNR {snr:.1f}", fontsize=9)

fig.supxlabel("k cross-track (cycles/km); white circles = 1000 m and 500 m wavelength")
fig.supylabel("k along-track (cycles/km)")
fig.suptitle(f"SWOT pass {PASS}: 2D spectra of {TILE_KM} km tiles (high-passed SSHA)")
plt.savefig(f"figures/swell_scan_{PASS}.png", dpi=130, bbox_inches="tight")
plt.show()
