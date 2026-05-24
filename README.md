# opt-ddc-wx

Tools for estimating Double Disc Court playable weather hours from ERA5 hourly
weather data.

## CONUS Playable-Hours Map

Download monthly ERA5 hourly grids for the conterminous US:

```sh
python3 scripts/download_era5_conus_by_month.py --year 2022 --output-dir data/output/conus
```

This requires the Python `cdsapi` package and a configured CDS API key.

Then compute and plot the number of daylight hours where temperature is above 40 F
and 10 m wind speed is below 5 m/s. The plot includes county, state, and
country borders and requires the R `maps` package.

```sh
Rscript scripts/plot_conus_optimal_hours.R \
  --input_dir data/output/conus \
  --year 2020 \
  --output_csv data/output/conus/conus_optimal_hours_2020.csv \
  --output_png data/output/conus/conus_optimal_hours_2020.png
```

The downloader accepts `--months`, `--area NORTH WEST SOUTH EAST`, and optional
`--grid LAT_RES LON_RES` arguments. The plotter accepts `--min_temp_f` and
`--max_wind_mps` if you want to adjust the playable-weather definition. Use
`--year` to plot one year when multiple years of monthly grids are in the input
directory.

## Point-Location Scripts

The older single-location scripts are kept in `scripts/` for reference and
comparison:

- `download_era5_hourly.py`
- `download_era5_hourly_by_month.py`
- `Optimal-DDC-wx.R`
- `Optimal-DDC-wx-by-hour.R`
