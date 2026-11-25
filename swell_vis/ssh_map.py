"""Map SWOT SSHA and absolute dynamic topography (L2 LR SSH Expert, 2 km) over the Gulf Stream."""
import earthaccess
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

BBOX = (-80, 30, -50, 45)  # lon_min, lat_min, lon_max, lat_max
DATES = ("2026-10-31", "2026-11-20")

earthaccess.login()
results = earthaccess.search_data(
    short_name="SWOT_L2_LR_SSH_D",
    granule_name="*Expert*",  # Expert (not Basic) includes mean_dynamic_topography
    bounding_box=BBOX,
    temporal=DATES,
)
print(f"Found {len(results)} granules")
files = earthaccess.download(results, "data/")

fig, axes = plt.subplots(2, 1, figsize=(11, 12), subplot_kw={"projection": ccrs.PlateCarree()})
for ax in axes:
    ax.set_extent([BBOX[0], BBOX[2], BBOX[1], BBOX[3]])
    ax.add_feature(cfeature.LAND, color="0.8", zorder=2)
    ax.coastlines(zorder=3)
    ax.gridlines(draw_labels=True, linewidth=0.3)

for f in files:
    ds = xr.open_dataset(f)
    # height_cor_xover is NOT included in ssha_karin_2: it removes KaRIn roll/phase
    # errors (~±1 m cross-track tilts) estimated from crossover calibration
    ssha = ds.ssha_karin_2 + ds.height_cor_xover
    ssha = ssha.where((ds.ssha_karin_2_qual == 0) & (ds.height_cor_xover_qual == 0))
    # ADT = SSHA + MDT: put the time-mean circulation (incl. the mean Gulf Stream) back in
    adt = ssha + ds.mean_dynamic_topography
    lon = (ds.longitude + 180) % 360 - 180
    # keep only the part of the pass inside the box (passes are half-orbit long)
    inside = ((ds.latitude > BBOX[1] - 1) & (ds.latitude < BBOX[3] + 1)).any("num_pixels")
    lon, lat = lon[inside], ds.latitude[inside]
    im0 = axes[0].pcolormesh(lon, lat, ssha[inside], cmap="RdBu_r", vmin=-0.5, vmax=0.5,
                             transform=ccrs.PlateCarree(), shading="auto")
    im1 = axes[1].pcolormesh(lon, lat, adt[inside], cmap="Spectral_r", vmin=-0.6, vmax=1.2,
                             transform=ccrs.PlateCarree(), shading="auto")

plt.colorbar(im0, ax=axes[0], label="SSHA (m)", shrink=0.8)
plt.colorbar(im1, ax=axes[1], label="ADT (m)", shrink=0.8)
axes[0].set_title(f"SSHA (anomaly from mean sea surface), {DATES[0]} to {DATES[1]}")
axes[1].set_title("ADT = SSHA + MDT (height above geoid)")
plt.savefig("figures/gulf_stream_ssha_adt.png", dpi=150, bbox_inches="tight")
plt.show()
