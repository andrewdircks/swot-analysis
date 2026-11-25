"""Gulf Stream ADT: gap-free L4 map (DUACS, all altimeters) with recent SWOT KaRIn passes on top."""
import copernicusmarine
import earthaccess
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

BBOX = (-80, 30, -50, 45)  # lon_min, lat_min, lon_max, lat_max
DAY = "2026-11-19"                      # L4 map date
SWOT_DATES = ("2026-11-17", "2026-11-19")  # SWOT passes close in time to the map

# L4: daily gridded ADT from all altimeters (optimal interpolation; ~1/8 deg grid, ~100 km effective resolution)
# Login once beforehand with `copernicusmarine login` (stores credentials in ~/.copernicusmarine)
l4 = copernicusmarine.open_dataset(
    dataset_id="cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D",
    variables=["adt"],
    minimum_longitude=BBOX[0], maximum_longitude=BBOX[2],
    minimum_latitude=BBOX[1], maximum_latitude=BBOX[3],
    start_datetime=DAY, end_datetime=DAY,
).adt.isel(time=0)

# SWOT L2 Expert passes (same processing as ssh_map.py)
earthaccess.login()
results = earthaccess.search_data(short_name="SWOT_L2_LR_SSH_D", granule_name="*Expert*",
                                  bounding_box=BBOX, temporal=SWOT_DATES)
print(f"Found {len(results)} SWOT granules")
files = earthaccess.download(results, "data/")

fig = plt.figure(figsize=(13, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([BBOX[0], BBOX[2], BBOX[1], BBOX[3]])
ax.add_feature(cfeature.LAND, color="0.8", zorder=2)
ax.coastlines(zorder=3)
ax.gridlines(draw_labels=True, linewidth=0.3)

style = dict(cmap="Spectral_r", vmin=-0.6, vmax=1.2, transform=ccrs.PlateCarree())
ax.pcolormesh(l4.longitude, l4.latitude, l4, alpha=0.6, **style)  # faded background

for f in files:
    ds = xr.open_dataset(f)
    adt = ds.ssha_karin_2 + ds.height_cor_xover + ds.mean_dynamic_topography
    adt = adt.where((ds.ssha_karin_2_qual == 0) & (ds.height_cor_xover_qual == 0))
    lon = (ds.longitude + 180) % 360 - 180
    inside = ((ds.latitude > BBOX[1] - 1) & (ds.latitude < BBOX[3] + 1)).any("num_pixels")
    im = ax.pcolormesh(lon[inside], ds.latitude[inside], adt[inside], shading="auto", **style)

# L4 contours = smoothed geostrophic streamlines; compare with the sharper SWOT swaths
ax.contour(l4.longitude, l4.latitude, l4, levels=[-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0],
           colors="k", linewidths=0.6, transform=ccrs.PlateCarree())

plt.colorbar(im, ax=ax, label="ADT (m)", shrink=0.7)
ax.set_title(f"ADT: L4 DUACS {DAY} (faded, contours every 0.2 m)\n"
             f"+ SWOT KaRIn passes {SWOT_DATES[0]} to {SWOT_DATES[1]} (vivid)")
plt.savefig("figures/gulf_stream_l4_swot.png", dpi=150, bbox_inches="tight")
plt.show()
