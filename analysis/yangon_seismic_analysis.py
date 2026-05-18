#!/usr/bin/env python3
"""
Yangon Mw5.2 Seismic Microzonation — Complete Analysis Pipeline

Station data source: Geographical Society of Myanmar
  Yangon Mw5.2 Shake Feel Level 1-5 Microzonation (KML dataset)

Analyses:
  1. Geospatial Map  (IDW surface + building footprint overlay)
  2. Ground Motion Attenuation  (PGA/PGV vs distance, site amplification)
  3. Response Spectrum Analysis  (per zone, amplification ratio)
  4. Building Risk Assessment  (risk matrix, building counts, priority ranking)
  5. Statistical Analysis  (correlation, zone stats, regression, outliers)
"""

import os
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import BoundaryNorm, ListedColormap
import seaborn as sns
from scipy import stats
from scipy.spatial import cKDTree
import geopandas as gpd
import contextily as ctx

warnings.filterwarnings("ignore")
np.random.seed(42)

# ── Paths (relative to project root, one level above this script) ─────────────
ROOT      = Path(__file__).resolve().parent.parent
KML_PATH  = ROOT / "data" / "yangon_mw52_microzonation.kml"
BLDG_PATH = ROOT / "yangon_townships_buildings.geojson"
FIG_DIR   = ROOT / "docs" / "assets" / "figures"
DATA_DIR  = ROOT / "data"
FIG_DIR.mkdir(parents=True, exist_ok=True)

LEVEL_COLORS = {1: "#3399ff", 2: "#66cc33", 3: "#ffcc00", 4: "#ff6600", 5: "#cc0000"}
ZONE_COLORS  = {
    "Zone 1 - Basin Core":       "#cc0000",
    "Zone 2 - Mixed Urban":      "#ff6600",
    "Zone 3 - Baseline/Stiffer": "#009900",
}
CRS = "EPSG:4326"


# ── 0. Parse KML ─────────────────────────────────────────────────────────────
def parse_kml(path):
    ns = {"k": "http://www.opengis.net/kml/2.2"}
    root = ET.parse(path).getroot()
    records = []
    for pm in root.findall(".//k:Placemark", ns):
        coords = pm.find(".//k:coordinates", ns).text.strip()
        lon, lat, _ = map(float, coords.split(","))
        rec = {"name": pm.find("k:name", ns).text, "lon": lon, "lat": lat}
        for d in pm.findall(".//k:Data", ns):
            key = d.get("name")
            val = d.find("k:value", ns).text
            try:
                val = float(val)
            except (ValueError, TypeError):
                pass
            rec[key] = val
        records.append(rec)
    df = pd.DataFrame(records)
    df["Shake_Feel_Level_1to5"] = df["Shake_Feel_Level_1to5"].astype(int)
    df["sp_id"] = df["name"].str.extract(r"(SP\d+)")
    return df


# ── Helper: IDW interpolation ─────────────────────────────────────────────────
def idw(x, y, z, xi, yi, power=2):
    shape = xi.shape
    xi_f, yi_f = xi.ravel(), yi.ravel()
    dx = xi_f[:, None] - x[None, :]
    dy = yi_f[:, None] - y[None, :]
    d  = np.sqrt(dx**2 + dy**2)
    d  = np.where(d == 0, 1e-12, d)
    w  = 1.0 / d**power
    return (np.sum(w * z[None, :], axis=1) / np.sum(w, axis=1)).reshape(shape)


# ── 0b. Load geodata ──────────────────────────────────────────────────────────
def load_geodata(df, bldg_path):
    gdf_stations = gpd.GeoDataFrame(
        df, geometry=gpd.points_from_xy(df.lon, df.lat), crs=CRS
    )

    if not os.path.exists(bldg_path):
        print(f"  Warning: {bldg_path} not found — skipping building overlay")
        return gdf_stations, None

    print("  >> Loading building footprints (filtering to study area)...")
    pad = 0.05
    bbox = (df.lon.min() - pad, df.lat.min() - pad,
            df.lon.max() + pad, df.lat.max() + pad)
    gdf_bldg = gpd.read_file(bldg_path, bbox=bbox)
    print(f"     {len(gdf_bldg):,} buildings in study area")

    if gdf_bldg.empty:
        return gdf_stations, None

    # Assign shake level by nearest station (KNN k=1)
    centroids = gdf_bldg.geometry.centroid
    bldg_pts  = np.c_[centroids.x, centroids.y]
    sta_pts   = np.c_[df.lon.values, df.lat.values]
    _, idx    = cKDTree(sta_pts).query(bldg_pts, k=1)
    gdf_bldg["shake_level"] = df.Shake_Feel_Level_1to5.iloc[idx].values

    return gdf_stations, gdf_bldg


# ── Analysis 1: Geospatial Map ────────────────────────────────────────────────
def analysis_geospatial(df, gdf_stations, gdf_bldg):
    print("  >> IDW surface + building overlay map")
    pad = 0.03
    lon_min = df.lon.min() - pad;  lon_max = df.lon.max() + pad
    lat_min = df.lat.min() - pad;  lat_max = df.lat.max() + pad

    xi = np.linspace(lon_min, lon_max, 300)
    yi = np.linspace(lat_min, lat_max, 300)
    XI, YI = np.meshgrid(xi, yi)
    ZI = idw(df.lon.values, df.lat.values,
             df.Shake_Feel_Level_1to5.values.astype(float), XI, YI)

    cmap   = ListedColormap(["#3399ff", "#66cc33", "#ffcc00", "#ff6600", "#cc0000"])
    bounds = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    norm   = BoundaryNorm(bounds, cmap.N)

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # --- Map A: IDW shake surface ---
    ax = axes[0]
    ax.contourf(XI, YI, ZI, levels=bounds, cmap=cmap, norm=norm, alpha=0.55)
    ax.contour(XI, YI, ZI, levels=bounds, colors="k", linewidths=0.4, alpha=0.4)
    try:
        ctx.add_basemap(ax, crs=CRS,
                        source=ctx.providers.OpenStreetMap.Mapnik, alpha=0.45, zoom=11)
    except Exception:
        pass
    for level in sorted(df.Shake_Feel_Level_1to5.unique(), reverse=True):
        sub = df[df.Shake_Feel_Level_1to5 == level]
        ax.scatter(sub.lon, sub.lat, c=LEVEL_COLORS[level], s=180,
                   edgecolors="black", linewidths=0.8, zorder=6, label=f"Level {level}")
        for _, r in sub.iterrows():
            ax.annotate(r.sp_id, (r.lon, r.lat), textcoords="offset points",
                        xytext=(4, 4), fontsize=6.5, fontweight="bold")
    ax.set_xlim(lon_min, lon_max);  ax.set_ylim(lat_min, lat_max)
    ax.set_xlabel("Longitude");     ax.set_ylabel("Latitude")
    ax.set_title("IDW Shake Intensity Surface + Stations", fontweight="bold")
    ax.legend(title="Level", loc="upper left", fontsize=8)

    # --- Map B: Building footprints coloured by shake level ---
    ax2 = axes[1]
    if gdf_bldg is not None:
        for level in sorted(LEVEL_COLORS):
            sub_b = gdf_bldg[gdf_bldg.shake_level == level]
            if not sub_b.empty:
                sub_b.plot(ax=ax2, color=LEVEL_COLORS[level],
                           linewidth=0, alpha=0.55, label=f"Level {level}")
    else:
        ax2.text(0.5, 0.5, "Building footprints\nnot loaded",
                 ha="center", va="center", transform=ax2.transAxes, fontsize=12)
    gdf_stations.plot(ax=ax2,
                      color=[LEVEL_COLORS[l] for l in df.Shake_Feel_Level_1to5],
                      markersize=100, edgecolor="black", linewidth=0.8, zorder=5)
    for _, r in df.iterrows():
        ax2.annotate(r.sp_id, (r.lon, r.lat), textcoords="offset points",
                     xytext=(3, 3), fontsize=6.5, fontweight="bold")
    try:
        ctx.add_basemap(ax2, crs=CRS,
                        source=ctx.providers.OpenStreetMap.Mapnik, alpha=0.35, zoom=11)
    except Exception:
        pass
    ax2.set_xlim(lon_min, lon_max);  ax2.set_ylim(lat_min, lat_max)
    ax2.set_xlabel("Longitude");      ax2.set_ylabel("Latitude")
    ax2.set_title("Building Footprints by Assigned Shake Level\n(nearest-station)", fontweight="bold")
    patches = [mpatches.Patch(color=LEVEL_COLORS[l], label=f"Level {l}")
               for l in sorted(LEVEL_COLORS)]
    ax2.legend(handles=patches, title="Level", loc="upper left", fontsize=8)

    plt.suptitle("Yangon Mw5.2 — Geospatial Shake Intensity Map",
                 fontweight="bold", fontsize=13)
    plt.tight_layout()
    out = str(FIG_DIR / "fig1_geospatial_map.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"     Saved: {out}")


# ── Analysis 2: Attenuation ───────────────────────────────────────────────────
def analysis_attenuation(df):
    print("  >> Attenuation analysis")
    zone_short = {k: k.split(" - ")[1] for k in ZONE_COLORS}
    r_range = np.linspace(df.Repi_km.min() * 0.9, df.Repi_km.max() * 1.05, 100)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    def scatter_zones(ax, xcol, ycol, xlabel, ylabel, title, logfit=True):
        for zone, color in ZONE_COLORS.items():
            sub = df[df.zone == zone]
            ax.scatter(sub[xcol], sub[ycol], c=color, s=80, edgecolors="k",
                       linewidths=0.5, label=zone_short[zone], zorder=4)
            if logfit and len(sub) > 2:
                try:
                    c = np.polyfit(np.log(sub[xcol]), np.log(sub[ycol]), 1)
                    ax.plot(r_range, np.exp(c[1]) * r_range ** c[0],
                            color=color, linestyle="--", linewidth=1.2, alpha=0.7)
                except Exception:
                    pass
        ax.set_xlabel(xlabel);  ax.set_ylabel(ylabel)
        ax.set_title(title, fontweight="bold")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    scatter_zones(axes[0, 0], "Repi_km", "M5p2_PGA_g",
                  "Epicentral Distance (km)", "PGA (g)",
                  "PGA Attenuation (dashed = power-law fit)")
    scatter_zones(axes[0, 1], "Repi_km", "M5p2_PGV_cm_s",
                  "Epicentral Distance (km)", "PGV (cm/s)",
                  "PGV Attenuation")

    # Site class amplification
    ax = axes[1, 0]
    site_pal = {"D - Soft Soil": "#3399ff", "E - Very Soft / Basin": "#cc0000"}
    for site, color in site_pal.items():
        sub = df[df.Hybrid_Site_Class == site]
        ax.scatter(sub.Repi_km, sub.M5p2_PGA_g, c=color, s=80, edgecolors="k",
                   linewidths=0.5, label=site, zorder=4)
    d_mean = df[df.Hybrid_Site_Class == "D - Soft Soil"].M5p2_PGA_g.mean()
    e_mean = df[df.Hybrid_Site_Class == "E - Very Soft / Basin"].M5p2_PGA_g.mean()
    ax.text(0.05, 0.95,
            f"Mean PGA  Site E / Site D = {e_mean/d_mean:.2f}×",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="lightyellow", ec="orange"))
    ax.set_xlabel("Epicentral Distance (km)");  ax.set_ylabel("PGA (g)")
    ax.set_title("Site Amplification: PGA by Site Class", fontweight="bold")
    ax.legend(fontsize=8);  ax.grid(True, alpha=0.3)

    # Shake level vs Repi (jittered)
    ax = axes[1, 1]
    jitter = np.random.uniform(-0.15, 0.15, len(df))
    for zone, color in ZONE_COLORS.items():
        sub = df[df.zone == zone]
        ax.scatter(sub.Repi_km,
                   sub.Shake_Feel_Level_1to5 + jitter[sub.index],
                   c=color, s=80, edgecolors="k", linewidths=0.5,
                   label=zone_short[zone], alpha=0.85)
    ax.set_xlabel("Epicentral Distance (km)");  ax.set_ylabel("Shake/Feel Level")
    ax.set_yticks([1, 2, 3, 4, 5])
    ax.set_title("Shake Level vs Distance (jittered)", fontweight="bold")
    ax.legend(fontsize=8);  ax.grid(True, alpha=0.3)

    plt.suptitle("Ground Motion Attenuation Analysis — Yangon Mw5.2",
                 fontweight="bold", fontsize=13)
    plt.tight_layout()
    out = str(FIG_DIR / "fig2_attenuation.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"     Saved: {out}")


# ── Analysis 3: Response Spectra ──────────────────────────────────────────────
def analysis_response_spectra(df):
    print("  >> Response spectrum analysis")
    periods = [0.3, 1.0, 2.0]
    sa_cols = ["M5p2_Sa_0p3_g", "M5p2_Sa_1p0_g", "M5p2_Sa_2p0_g"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # --- 3a: Individual spectra by zone ---
    ax = axes[0]
    for _, row in df.iterrows():
        ax.plot(periods, [row[c] for c in sa_cols],
                color=ZONE_COLORS.get(row.zone, "gray"), alpha=0.35, linewidth=1)
    for zone, color in ZONE_COLORS.items():
        sub = df[df.zone == zone]
        mean_sa = [sub[c].mean() for c in sa_cols]
        ax.plot(periods, mean_sa, color=color, linewidth=2.8, zorder=5,
                label=f"{zone.split(' - ')[1]} (mean)")
    ax.set_xlabel("Period (s)");  ax.set_ylabel("Sa (g)")
    ax.set_title("Response Spectra by Station\n(bold = zone mean)", fontweight="bold")
    ax.legend(fontsize=8);  ax.grid(True, alpha=0.3)

    # --- 3b: Boxplot Sa by period and zone ---
    ax = axes[1]
    rows = []
    for col, T in zip(sa_cols, periods):
        for _, r in df.iterrows():
            rows.append({"Period": f"{T}s", "Sa (g)": r[col],
                         "Zone": r.zone.split(" - ")[1]})
    pdf = pd.DataFrame(rows)
    pal = {"Basin Core": "#cc0000", "Mixed Urban": "#ff6600",
           "Baseline/Stiffer": "#009900"}
    sns.boxplot(data=pdf, x="Period", y="Sa (g)", hue="Zone", palette=pal, ax=ax)
    ax.set_title("Sa Distribution by Period & Zone", fontweight="bold")
    ax.legend(fontsize=8);  ax.grid(True, alpha=0.3, axis="y")

    # --- 3c: Amplification ratio Zone 1 / Zone 3 ---
    ax = axes[2]
    z1 = df[df.zone == "Zone 1 - Basin Core"][sa_cols].mean()
    z3 = df[df.zone == "Zone 3 - Baseline/Stiffer"][sa_cols].mean()
    ratios = (z1 / z3).values
    bars = ax.bar([f"{T}s" for T in periods], ratios,
                  color=["#cc3300", "#cc6600", "#cc9900"], edgecolor="black")
    ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8)
    for bar, val in zip(bars, ratios):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.05,
                f"{val:.2f}×", ha="center", fontsize=11, fontweight="bold")
    ax.set_ylabel("Amplification Ratio (Zone 1 / Zone 3)")
    ax.set_title("Basin Amplification Factor\n(Zone 1 vs Zone 3)", fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")

    plt.suptitle("Response Spectrum Analysis — Yangon Mw5.2",
                 fontweight="bold", fontsize=13)
    plt.tight_layout()
    out = str(FIG_DIR / "fig3_response_spectra.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"     Saved: {out}")


# ── Analysis 4: Building Risk ─────────────────────────────────────────────────
def analysis_building_risk(df, gdf_bldg):
    print("  >> Building risk assessment")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # --- 4a: Risk matrix Zone × Level ---
    ax = axes[0]
    rm = df.groupby(["zone", "Shake_Feel_Level_1to5"]).size().unstack(fill_value=0)
    zlabels = [z.split(" - ")[1] for z in rm.index]
    im = ax.imshow(rm.values, cmap="Reds", aspect="auto", vmin=0)
    ax.set_xticks(range(len(rm.columns)))
    ax.set_xticklabels([f"L{c}" for c in rm.columns])
    ax.set_yticks(range(len(zlabels)))
    ax.set_yticklabels(zlabels)
    ax.set_title("Risk Matrix: Zone × Shake Level\n(station count)", fontweight="bold")
    for i in range(rm.shape[0]):
        for j in range(rm.shape[1]):
            v = rm.values[i, j]
            if v > 0:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=14,
                        fontweight="bold", color="white" if v > 1 else "black")
    plt.colorbar(im, ax=ax, shrink=0.7)

    # --- 4b: Building count by shake level ---
    ax = axes[1]
    if gdf_bldg is not None:
        bldg_counts = gdf_bldg.groupby("shake_level").size()
        all_levels  = pd.Series(0, index=range(1, 6))
        all_levels.update(bldg_counts)
        colors = [LEVEL_COLORS[l] for l in all_levels.index]
        bars = ax.bar(all_levels.index, all_levels.values,
                      color=colors, edgecolor="black")
        for bar, val in zip(bars, all_levels.values):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(all_levels.values) * 0.01,
                    f"{val:,}", ha="center", fontsize=9, fontweight="bold")
        ax.set_xlabel("Shake/Feel Level"); ax.set_ylabel("Building Count")
        ax.set_title("Buildings by Assigned Shake Level\n(nearest-station)", fontweight="bold")
        ax.set_xticks(range(1, 6)); ax.grid(True, alpha=0.3, axis="y")
    else:
        ax.text(0.5, 0.5, "Building footprints\nnot available",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_title("Buildings by Shake Level", fontweight="bold")

    # --- 4c: Priority ranking ---
    ax = axes[2]
    df_s   = df.sort_values("Shake_Feel_Level_1to5", ascending=True)
    sp_ids = df_s.sp_id.values
    levels = df_s.Shake_Feel_Level_1to5.values
    colors = [LEVEL_COLORS[l] for l in levels]
    ax.barh(range(len(sp_ids)), levels, color=colors, edgecolor="black", linewidth=0.4)
    ax.set_yticks(range(len(sp_ids))); ax.set_yticklabels(sp_ids, fontsize=8)
    ax.axvline(3.5, color="orange", linestyle="--", linewidth=1.2,
               label="Action threshold (Level > 3)")
    ax.set_xlabel("Shake/Feel Level"); ax.set_xlim(0, 6)
    ax.set_title("Station Priority Ranking", fontweight="bold")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3, axis="x")

    plt.suptitle("Building Risk Assessment — Yangon Mw5.2",
                 fontweight="bold", fontsize=13)
    plt.tight_layout()
    out = str(FIG_DIR / "fig4_building_risk.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"     Saved: {out}")


# ── Analysis 5: Statistical Analysis ──────────────────────────────────────────
def analysis_statistics(df):
    print("  >> Statistical analysis")
    num_cols = ["Shake_Feel_Level_1to5", "M5p2_PGA_g", "M5p2_PGV_cm_s",
                "M5p2_Sa_0p3_g", "M5p2_Sa_1p0_g", "M5p2_Sa_2p0_g", "Repi_km"]
    labels   = ["Shake\nLevel", "PGA", "PGV", "Sa\n0.3s", "Sa\n1.0s", "Sa\n2.0s", "Repi"]

    fig = plt.figure(figsize=(18, 12))
    gs  = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.35)

    # --- 5a: Correlation heatmap ---
    ax1 = fig.add_subplot(gs[0, 0:2])
    corr = df[num_cols].corr()
    mask = np.triu(np.ones_like(corr, dtype=bool))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="RdYlGn", center=0,
                ax=ax1, xticklabels=labels, yticklabels=labels,
                linewidths=0.5, vmin=-1, vmax=1, mask=mask)
    ax1.set_title("Correlation Matrix (lower triangle)", fontweight="bold")

    # --- 5b: Zone statistics table ---
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.axis("off")
    zs = df.groupby("zone").agg(
        N=("name", "count"),
        Avg_Level=("Shake_Feel_Level_1to5", "mean"),
        Avg_PGA=("M5p2_PGA_g", "mean"),
        Avg_PGV=("M5p2_PGV_cm_s", "mean"),
        Avg_Repi=("Repi_km", "mean"),
    ).round(3)
    zs.index = [z.split(" - ")[1] for z in zs.index]
    tbl = ax2.table(
        cellText=zs.values, rowLabels=zs.index,
        colLabels=["N", "Avg Lvl", "Avg PGA", "Avg PGV", "Avg Repi"],
        cellLoc="center", loc="center", bbox=[0, 0, 1, 1],
    )
    tbl.auto_set_font_size(False); tbl.set_fontsize(8)
    for j in range(len(zs.columns)):
        tbl[(0, j)].set_facecolor("#cccccc")
    ax2.set_title("Zone Statistics Summary", fontweight="bold", pad=30)

    # --- 5c: PGA violin by zone ---
    ax3 = fig.add_subplot(gs[1, 0])
    df2 = df.copy()
    df2["Zone"] = df2["zone"].str.split(" - ").str[1]
    pal = {"Basin Core": "#cc0000", "Mixed Urban": "#ff6600",
           "Baseline/Stiffer": "#009900"}
    sns.violinplot(data=df2, x="Zone", y="M5p2_PGA_g", palette=pal,
                   ax=ax3, inner="point")
    ax3.set_xlabel("");  ax3.set_ylabel("PGA (g)")
    ax3.set_title("PGA Distribution by Zone", fontweight="bold")
    ax3.grid(True, alpha=0.3, axis="y")

    # --- 5d: Repi vs Shake Level with regression ---
    ax4 = fig.add_subplot(gs[1, 1])
    sl, ic, r, p, _ = stats.linregress(df.Repi_km, df.Shake_Feel_Level_1to5)
    x_r = np.linspace(df.Repi_km.min(), df.Repi_km.max(), 100)
    ax4.scatter(df.Repi_km, df.Shake_Feel_Level_1to5,
                c=[LEVEL_COLORS[l] for l in df.Shake_Feel_Level_1to5],
                s=80, edgecolors="k", linewidths=0.5)
    ax4.plot(x_r, sl * x_r + ic, "k--", linewidth=1.5, label=f"r={r:.2f}, p={p:.3f}")
    ax4.set_xlabel("Epicentral Distance (km)"); ax4.set_ylabel("Shake/Feel Level")
    ax4.set_yticks([1, 2, 3, 4, 5])
    ax4.set_title("Repi vs Shake Level\n(Does distance explain shaking?)", fontweight="bold")
    significance = "Significant" if p < 0.05 else "Not significant"
    ax4.text(0.05, 0.05, f"r = {r:.2f}  p = {p:.3f}\n{significance} at α = 0.05",
             transform=ax4.transAxes, va="bottom", fontsize=9,
             bbox=dict(boxstyle="round", fc="lightyellow", ec="orange"))
    ax4.legend(fontsize=8);  ax4.grid(True, alpha=0.3)

    # --- 5e: Outlier detection (z-score on PGA + Repi) ---
    ax5 = fig.add_subplot(gs[1, 2])
    df3 = df.copy()
    df3["z_PGA"]  = stats.zscore(df3.M5p2_PGA_g)
    df3["z_Repi"] = stats.zscore(df3.Repi_km)
    df3["outlier"] = (df3.z_PGA.abs() >= 2) | (df3.z_Repi.abs() >= 2)

    normal   = df3[~df3.outlier]
    outliers = df3[df3.outlier]
    ax5.scatter(normal.Repi_km, normal.M5p2_PGA_g, c="steelblue", s=80,
                edgecolors="k", linewidths=0.5, label="Normal", zorder=3)
    if not outliers.empty:
        ax5.scatter(outliers.Repi_km, outliers.M5p2_PGA_g, c="red", s=150,
                    edgecolors="k", linewidths=0.8, marker="*",
                    label=f"Outlier |z|≥2 (n={len(outliers)})", zorder=4)
        for _, row in outliers.iterrows():
            ax5.annotate(row.sp_id, (row.Repi_km, row.M5p2_PGA_g),
                         textcoords="offset points", xytext=(5, 5), fontsize=8)
    ax5.set_xlabel("Repi (km)");  ax5.set_ylabel("PGA (g)")
    ax5.set_title("Outlier Detection\n(z-score on PGA and Repi)", fontweight="bold")
    ax5.legend(fontsize=8);  ax5.grid(True, alpha=0.3)

    plt.suptitle("Statistical Analysis — Yangon Mw5.2 Microzonation",
                 fontweight="bold", fontsize=13)
    out = str(FIG_DIR / "fig5_statistics.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"     Saved: {out}")


# ── Export CSV ────────────────────────────────────────────────────────────────
def export_csv(df):
    cols = [
        "sp_id", "name", "lon", "lat", "zone", "Hybrid_Site_Class",
        "Shake_Feel_Level_1to5", "Hybrid_Resonance_Hazard", "Today_MMI_est",
        "M5p2_PGA_g", "M5p2_PGV_cm_s",
        "M5p2_Sa_0p3_g", "M5p2_Sa_1p0_g", "M5p2_Sa_2p0_g",
        "Repi_km", "MNBC2025_Application_Zone", "Risk",
    ]
    out = str(DATA_DIR / "stations_data.csv")
    df[cols].to_csv(out, index=False)
    print(f"     Saved: {out}")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 58)
    print("  Yangon Mw5.2 Seismic Microzonation — Analysis Pipeline")
    print("=" * 58)

    print("\n[0] Parsing KML...")
    df = parse_kml(KML_PATH)
    print(f"    {len(df)} stations loaded")

    print("\n[0b] Loading geodata...")
    gdf_stations, gdf_bldg = load_geodata(df, BLDG_PATH)

    print("\n[1] Geospatial Map...")
    analysis_geospatial(df, gdf_stations, gdf_bldg)

    print("\n[2] Attenuation Analysis...")
    analysis_attenuation(df)

    print("\n[3] Response Spectrum Analysis...")
    analysis_response_spectra(df)

    print("\n[4] Building Risk Assessment...")
    analysis_building_risk(df, gdf_bldg)

    print("\n[5] Statistical Analysis...")
    analysis_statistics(df)

    print("\n[6] Exporting CSV...")
    export_csv(df)

    print(f"\nDone! All outputs saved in: {os.path.abspath(OUTPUT_DIR)}")
