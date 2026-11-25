"""Find the big swell: map SWOT KaRIn significant wave height (L2 WindWave, 2 km) in the NE Pacific, one panel per day."""
import earthaccess
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
import cartopy.feature as cfeature

BBOX = (-160, 25, -120, 55)  # lon_min, lat_min, lon_max, lat_max
DATES = ("2026-11-17", "2026-11-23")
BUOY = (-123.47, 37.94)  # Point Reyes, NDBC 46214

earthaccess.login()
results = earthaccess.search_data(short_name="SWOT_L2_LR_SSH_D", granule_name="*WindWave*",
                                  bounding_box=BBOX, temporal=DATES)
print(f"Found {len(results)} granules")
files = sorted(earthaccess.download(results, "data/"))

# group passes by day: file names look like SWOT_L2_LR_SSH_WindWave_059_164_20261118T045413_...
names = [str(f).split("WindWave_")[1] for f in files]
days = sorted({n[8:16] for n in names})
ncols = (len(days) + 1) // 2
fig, axes = plt.subplots(2, ncols, figsize=(5.5 * ncols, 9), layout="constrained",
                         subplot_kw={"projection": ccrs.PlateCarree()})
axes = axes.ravel()

for ax, day in zip(axes, days):
    ax.set_extent([BBOX[0], BBOX[2], BBOX[1], BBOX[3]])
    ax.add_feature(cfeature.LAND, color="0.8", zorder=2)
    ax.coastlines(zorder=3)
    ax.gridlines(draw_labels=["left", "bottom"], linewidth=0.3)
    ax.plot(*BUOY, "c*", markersize=14, transform=ccrs.PlateCarree(), zorder=5)
    ax.set_title(f"{day[6:8]} Nov {day[:4]}")

    for f, name in zip(files, names):
        if name[8:16] != day:
            continue
        ds = xr.open_dataset(f)
        # qual is a bit field: values < 2**17 are only "suspect" bits; drop "degraded" and "bad".
        # also drop the coastal zone, where land contaminates the radar return (gives bogus ~30 m SWH)
        swh = ds.swh_karin.where((ds.swh_karin_qual < 2**17) & (ds.distance_to_coast > 50e3))
        lon = (ds.longitude + 180) % 360 - 180
        inside = ((ds.latitude > BBOX[1]) & (ds.latitude < BBOX[3]) & (lon > BBOX[0]) & (lon < BBOX[2])).any("num_pixels")
        lon, lat, swh = lon[inside], ds.latitude[inside], swh[inside]
        if swh.isnull().all():
            continue
        im = ax.pcolormesh(lon, lat, swh, cmap="magma_r", vmin=0, vmax=10,
                           transform=ccrs.PlateCarree(), shading="auto")

        # label near 51N on the pass's ocean track, where crossing passes (an "X") are far apart: "pass, hour UTC"
        nlon, nlat = lon[:, 34].values, lat[:, 34].values  # nadir track (middle pixel)
        ocean = np.flatnonzero(swh.notnull().any("num_pixels").values)
        j = ocean[np.argmin(abs(nlat[ocean] - 51))]
        ax.text(nlon[j], nlat[j], f"{name[4:7]}\n{name[17:19]}Z", fontsize=9, ha="center", va="center",
                bbox=dict(boxstyle="round", fc="white", ec="k", lw=0.5), transform=ccrs.PlateCarree(), zorder=6)

        # peak SWH: median across each line, smoothed along-track over ~20 km
        rows = swh.median("num_pixels").rolling(num_lines=10, center=True).median()
        i = int(np.nanargmax(rows.values))
        halo = [pe.withStroke(linewidth=3, foreground="white")]  # white outline: readable on any colour
        ax.plot(nlon[i], nlat[i], "kx", ms=9, mew=2, path_effects=halo, transform=ccrs.PlateCarree(), zorder=5)
        ax.text(nlon[i] + 0.7, nlat[i], f"{float(rows[i]):.1f} m", fontsize=9, fontweight="bold",
                path_effects=halo, transform=ccrs.PlateCarree(), zorder=5)
        print(f"{name[:30]}  peak SWH {float(rows[i]):4.1f} m at {nlat[i]:5.1f}N {nlon[i]:6.1f}E")

for ax in axes[len(days):]:
    ax.remove()
fig.colorbar(im, ax=axes[:len(days)], label="SWH (m)", shrink=0.6)
fig.suptitle("SWOT KaRIn significant wave height. Labels: pass number + time (UTC); "
             "× = peak SWH on that pass; ★ = Point Reyes buoy 46214", fontsize=12)
plt.savefig("figures/swell_finder.png", dpi=120, bbox_inches="tight")
plt.show()
