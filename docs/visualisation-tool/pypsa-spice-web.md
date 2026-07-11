# PyPSA-SPICE power output explorer

The output explorer is the first web-based replacement slice for the Streamlit visualisation interface. It reads the existing result CSV structure directly and does not modify model outputs. The current release covers Power outputs and scenario comparison. Comparison mode places both scenarios side by side for each indicator; the optional difference view calculates `comparison − primary` for every aligned series and timestamp or model year.

## Run locally

From the repository root:

```bash
python -m pip install -r requirements-web.txt
npm install --prefix pypsa_spice_web/frontend
./run_web.sh
```

Open <http://127.0.0.1:8000>. The launcher uses the active Python environment when the web dependencies are installed, and falls back to the existing `hotpot` Conda environment when available.

Set `HOST` or `PORT` to override the defaults:

```bash
HOST=0.0.0.0 PORT=8080 ./run_web.sh
```

## Result layout

The app discovers data in the established structure:

```text
data/<dataset>/<project>/results/<scenario>/csvs/<sector>/<year>/<table>.csv
```

Yearly visualisations combine the tables found across modelled-year folders. Hourly visualisations use the selected year. Large hourly results are sampled for display while CSV downloads always return the full source file.

## Current scope

This release includes Power output visualisation and scenario comparison only. Model configuration, execution, Industry, Transport, Emissions, and Costs remain intentionally unavailable until their backend workflows are implemented.

The current data-management strategy is local discovery of the checked-out `data/` tree. Uploads, object storage, authentication, and hosted data persistence are deliberately deferred until a deployment target is selected.
