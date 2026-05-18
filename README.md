# Yangon Mw 5.2 Seismic Microzonation

Preliminary microzonation and shake/feel level analysis for the Greater Yangon
Metropolitan Area following an **Mw 5.2 earthquake south of Yangon, Myanmar**.

**Live report →** https://geonet-myanmar.github.io/yangon-earthquake-analysis/

---

## What this project does

Classifies 22 monitoring stations into five **Shake/Feel Levels (L1–L5)** by
integrating:

- Peak Ground Acceleration (PGA) and Velocity (PGV) from GMPE attenuation
- Spectral accelerations at T = 0.3, 1.0, and 2.0 s
- Hybrid site-class assignment (V<sub>s30</sub>-based: Site Class D / E)
- Resonance hazard assessment for the Yangon alluvial basin
- Cross-reference against **MNBC 2025** application zones

Outputs include IDW-interpolated shake maps, attenuation curves, response spectra,
building risk counts (426,661 footprints), and a full statistical analysis.

---

## Repository structure

```
.
├── analysis/
│   ├── yangon_seismic_analysis.py   # Complete Python analysis pipeline
│   └── requirements.txt             # Python dependencies
├── data/
│   ├── yangon_mw52_microzonation.kml  # Source station data — Geographical Society of Myanmar
│   └── stations_data.csv              # Processed station table
├── docs/                            # GitHub Pages site root
│   ├── index.html                   # Report website (all sections)
│   ├── css/style.css
│   ├── js/main.js
│   └── assets/figures/              # Generated analysis figures (PNG)
│       ├── fig1_geospatial_map.png
│       ├── fig2_attenuation.png
│       ├── fig3_response_spectra.png
│       ├── fig4_building_risk.png
│       └── fig5_statistics.png
└── results/                         # Raw analysis output (gitignored)
```

> **Note:** `yangon_townships_buildings.geojson` (236 MB) is excluded from the
> repository by `.gitignore`. Obtain the building footprint dataset separately
> and place it in the project root before running the analysis script.

---

## Quick start

### 1. Install dependencies

```bash
pip install -r analysis/requirements.txt
```

### 2. Run the full analysis pipeline

```bash
# From the project root
python analysis/yangon_seismic_analysis.py
```

This reads from `data/`, computes all analyses, and writes figures to
`docs/assets/figures/` and a CSV to `data/`.

**Building footprints (optional):** Place `yangon_townships_buildings.geojson`
in the project root to enable the building overlay and count analysis.

### 3. Preview the report locally

Open `docs/index.html` in any modern browser — no web server required.

---

## Analysis components

| Figure | Description |
|--------|-------------|
| **Fig 1** | IDW-interpolated shake intensity surface + 426 K building footprint overlay |
| **Fig 2** | PGA/PGV attenuation vs distance; site class amplification comparison |
| **Fig 3** | Response spectra by zone; basin amplification factor by period |
| **Fig 4** | Risk matrix (Zone × Level); building counts; station priority ranking |
| **Fig 5** | Correlation matrix; zone stats; regression (Repi vs Level); outlier detection |

---

## Key results

| Zone | Site Class | Shake Levels | MNBC 2025 |
|------|-----------|--------------|-----------|
| Zone 1 – Basin Core | E (Very Soft, V<sub>s30</sub> < 180 m/s) | **L4 – L5** | Zone 3 (site-specific spectrum) |
| Zone 2 – Mixed Urban | D (Soft Soil) | L3 – L4 | Zone 2 |
| Zone 3 – Baseline/Stiffer | D (Soft Soil) | L2 – L3 | Zone 2 |

Basin amplification factor (Zone 1 / Zone 3):
- **T = 0.3 s → 1.68×**
- **T = 1.0 s → 5.75×**
- **T = 2.0 s → 7.97×**

---

## Deploying to GitHub Pages

The site is deployed at:
**https://geonet-myanmar.github.io/yangon-earthquake-analysis/**

Source: `main` branch, `/docs` folder (GitHub Pages setting).

---

## Dependencies

See `analysis/requirements.txt`. Core packages:

```
pandas >= 2.0
numpy >= 1.24
matplotlib >= 3.7
seaborn >= 0.12
scipy >= 1.10
geopandas >= 0.13
contextily >= 1.3
shapely >= 2.0
```

Python 3.10 or later recommended.

---

## Data attribution

Station data (KML) provided by the **Geographical Society of Myanmar**.
Please credit the Geographical Society of Myanmar in any publication or
derivative work that uses this dataset.

---

## Disclaimer

This is a **preliminary simulation** based on GMPE predictions and proxy-based
site classification. Results should be validated with instrumental recordings
and geotechnical borehole data before use in formal structural design. The
authors assume no liability for decisions made solely on the basis of this report.
