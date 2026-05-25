---

editor_options: 
  markdown: 
    wrap: 72
---

# opt-ddc-wx

Tools for estimating Double Disc Court playable weather hours from ERA5 hourly weather data.

## CONUS Playable-Hours Map

Download monthly ERA5 hourly grids for the conterminous US:

``` sh
python3 scripts/download_era5_conus_by_month.py --year 2025 --output-dir data/output/conus
```

This requires the Python `cdsapi` package and a configured CDS API key.

Then compute and plot the number of daylight hours where temperature is above 40 F and 10 m wind speed is below 5 m/s. The plot includes county, state, and country borders and requires the R `maps` package.

``` sh
Rscript scripts/plot_conus_optimal_hours.R \
  --input_dir data/output/conus \
  --year 2025 \
  --output_csv data/output/conus/conus_optimal_hours_2025.csv \
  --output_png data/output/conus/conus_optimal_hours_2025.png
```

The downloader accepts `--months`, `--area NORTH WEST SOUTH EAST`, and optional `--grid LAT_RES LON_RES` arguments. The plotter accepts `--min_temp_f` and `--max_wind_mps` if you want to adjust the playable-weather definition. Use `--year` to plot one year when multiple years of monthly grids are in the input directory.

After generating annual summary CSVs, average multiple years and plot only the lower and upper 10% of CONUS grid cells:

``` sh
Rscript scripts/plot_conus_optimal_percentiles.R \
  --input_dir data/output/conus \
  --years 2020:2025 \
  --output_csv data/output/conus/conus_optimal_hours_2020_2025_average.csv \
  --output_png data/output/conus/conus_optimal_hours_2020_2025_percentiles.png
```

To repeat the percentile classification separately for every calendar month and
make a 12-panel map of the lower and upper 20% of CONUS grid cells:

```sh
Rscript scripts/plot_conus_monthly_optimal_percentiles.R \
  --input_dir data/output/conus \
  --years 2020:2025 \
  --output_csv data/output/conus/conus_monthly_optimal_hours_2020_2025_average.csv \
  --output_png data/output/conus/conus_monthly_optimal_hours_2020_2025_percentiles.png
```

## Point-Location Scripts

The older single-location scripts are kept in `scripts/` for reference and comparison:

- `download_era5_hourly.py`
- `download_era5_hourly_by_month.py`
- `Optimal-DDC-wx.R`
- `Optimal-DDC-wx-by-hour.R`
