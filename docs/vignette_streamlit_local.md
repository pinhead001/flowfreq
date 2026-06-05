# Vignette: Run the Streamlit App Locally

The HydroLib Streamlit app (`app/streamlit_app.py`) provides an interactive browser-based interface for downloading USGS streamflow data, running Bulletin 17C analysis, and exporting results.

## 1. Install

```bash
# From the repo root
cd C:/a/hal/hybrid-17c-cld   # adjust to your clone path

python -m venv .venv
.venv\Scripts\activate        # Windows
# source .venv/bin/activate   # macOS/Linux

pip install -e ".[dev]"
pip install streamlit
```

Verify:
```bash
python -c "import streamlit, hydrolib; print('ready')"
```

## 2. Run

**Always run from the repo root** — the app uses relative imports (`from app.ffa_runner import ...`):

```bash
cd C:/a/hal/hybrid-17c-cld
streamlit run app/streamlit_app.py
```

The browser opens automatically at `http://localhost:8501`.

To specify a port:
```bash
streamlit run app/streamlit_app.py --server.port 8502
```

## 3. Using the App

### Daily Flow Mode (default)

1. **Input Mode** — choose Single Gage or Multiple Gages
2. Enter a USGS gage number (e.g., `03606500` for Big Sandy River, TN)
3. **Plot Options** — check Daily Time Series, Summary Hydrograph, Flow Duration Curve
4. Click **Download Data** — fetches from NWIS, generates plots
5. Use the **Plot Date Range** sliders to filter the plotted period without re-downloading
6. Click **Update Plots** to redraw

### Flood Frequency Analysis

1. Check **Enable Flood Frequency Analysis** in the sidebar
2. Set **Regional Skew** (default −0.302, nationwide B17C mean) and **Regional Skew SE** (default 0.55)
3. **Skew Options** — check any combination of:
   - **Station Skew** — raw skew estimated from the data
   - **Weighted Skew** *(default)* — B17C weighted combination of station + regional
   - **Regional Skew** — applies the regional value directly
4. Check **Frequency Curve** to display the LP3 plot
5. Click **Download Data** — downloads peak flows and runs EMA analysis
6. View results in the **Flood Frequency Results** expander below each gage:
   - Convergence badge (EMA Converged / MOM Fallback)
   - LP3 Parameters table (μ, σ, skews)
   - One frequency table per selected skew (RI 1.5–500 yr with 90% CI)

### Multi-Gage Mode

Enter multiple gage numbers (one per line), enable FFA, and a **Flood Frequency Comparison** table appears comparing 100-yr flows across all gages.

### Export

1. Sidebar → **Export** section → select gages
2. Click **Download ZIP** — file contains per-gage subdirectories:
   ```
   03606500/
     daily_flow.csv
     daily_timeseries.png
     summary_hydrograph.png
     flow_duration_curve.png
     frequency_curve.png         (if FFA enabled)
     frequency_table.csv         (9 return intervals)
     lp3_parameters.csv
   comparison_summary.csv        (multi-gage only)
   ```

## 4. Troubleshooting

| Symptom | Fix |
|---------|-----|
| `ModuleNotFoundError: hydrolib` | Run from repo root; ensure venv is active |
| `ModuleNotFoundError: app.ffa_runner` | Run `streamlit run app/streamlit_app.py` from repo root, not from inside `app/` |
| Port already in use | `streamlit run app/streamlit_app.py --server.port 8502` |
| Blank page / no data | Check NWIS connectivity; some gages have restricted peak data |
| EMA did not converge | Short records (<10 peaks) trigger MOM fallback automatically |
