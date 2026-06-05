# Vignette: Deploy to Streamlit Community Cloud

Streamlit Community Cloud provides free hosting for public Streamlit apps linked to a public GitHub repository. Deployment takes about 5 minutes.

## Prerequisites

- A public GitHub repository containing the hydrolib code
- A free account at [share.streamlit.io](https://share.streamlit.io)

## 1. Prepare the Repository

The app needs two files at the repo root that Streamlit Cloud reads automatically.

### `app/requirements.txt` (already present)

```
streamlit>=1.28.0
numpy>=1.20.0
pandas>=1.3.0
matplotlib>=3.4.0
scipy>=1.7.0
requests>=2.25.0
```

Streamlit Cloud installs from this file automatically. The `hydrolib` package itself is installed from the repo via `pip install -e .` if a `setup.py` or `pyproject.toml` is present — or you can add it explicitly:

### `requirements.txt` (repo root — Streamlit Cloud primary)

Create or update `requirements.txt` at the repo root:

```
numpy>=1.20.0
pandas>=1.3.0
matplotlib>=3.4.0
scipy>=1.7.0
requests>=2.25.0
click>=8.0
streamlit>=1.28.0
```

Streamlit Cloud will also run `pip install -e .` from the repo root if `setup.py` or `pyproject.toml` is present, installing `hydrolib` as a package.

### `packages.txt` (optional — system packages)

Only needed if you use Fortran extensions (`hydrolib.peakfqr`). The compiled `.pyd` / `.so` files included in the repo should work on Linux (Streamlit Cloud runs Ubuntu). If not, add build tools:

```
gfortran
```

## 2. Push to GitHub

```bash
git push origin main
```

## 3. Deploy on Streamlit Cloud

1. Go to [share.streamlit.io](https://share.streamlit.io) and sign in with GitHub
2. Click **New app**
3. Fill in:
   - **Repository:** `your-username/hydrolib`
   - **Branch:** `main`
   - **Main file path:** `app/streamlit_app.py`
4. Click **Deploy**

Streamlit Cloud will:
- Clone the repo
- Install dependencies from `requirements.txt`
- Run `pip install -e .` (installs `hydrolib`)
- Start the app at `https://your-username-hydrolib-app-streamlit-app-XXXXXX.streamlit.app`

## 4. App URL and Settings

After deployment, the URL is shown in the Streamlit Cloud dashboard. You can:
- Set a **custom subdomain** in app settings (e.g., `hydrolib.streamlit.app`)
- **Reboot** the app if it goes to sleep (free tier sleeps after inactivity)
- View **logs** in the cloud dashboard for debugging

## 5. Updating the App

Push to `main` — Streamlit Cloud auto-deploys on every push:

```bash
git push origin main
# App redeploys in ~60 seconds
```

## 6. Environment Variables / Secrets

The app makes no authenticated API calls (USGS NWIS is public). No secrets are needed. If you add authenticated services later, use Streamlit Cloud's **Secrets** tab:

```toml
# .streamlit/secrets.toml  (local only — do NOT commit)
[api_keys]
some_key = "abc123"
```

Access in code:
```python
import streamlit as st
key = st.secrets["api_keys"]["some_key"]
```

## 7. Known Limitations on Streamlit Cloud

| Issue | Notes |
|-------|-------|
| Compiled Fortran (`.pyd`/`.so`) | `.pyd` files (Windows) won't run on Linux Cloud; `.so` files compiled on Linux will work. Add `gfortran` to `packages.txt` and a build step if needed. |
| Memory limits | Free tier: ~1 GB RAM. Large batch runs (many gages) may hit limits. |
| Sleep after inactivity | Free tier apps sleep after ~7 days of no traffic; wake on first visit (~30 s). |
| NWIS network access | Streamlit Cloud has outbound internet — USGS NWIS calls work fine. |
