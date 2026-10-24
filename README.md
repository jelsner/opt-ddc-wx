---

editor_options: 
  markdown: 
    wrap: 72
---

# opt-ddc-wx

Tools for estimating Double Disc Court playable weather hours from ERA5 hourly weather data.

## Tallahassee probability of wind at or below a threshold

This workflow uses station observations, independently of the ERA5 scripts.
Run from the project root; no R packages are required:

```r
source("scripts/tallahassee_wind_exceedance.R")

# A full table at 1 mph increments, for 09:00 through 15:00 local time.
p <- tlh_wind_exceedance("06-15")
head(p)
write.csv(p, "data/output/tlh_june15.csv", row.names = FALSE)

# Plot P(W > threshold) across local hours, with one line every 5 mph.
plot_tlh_wind(p)
plot_tlh_wind(p, file = "data/output/tlh_june15.png")
# Or calculate and plot together (still returns the probability table):
p <- tlh_wind_exceedance("06-15", plot = TRUE)

# Optional: plot P(wind <= threshold) instead.
plot_tlh_wind(p, type = "at_or_below")

# Display selected thresholds as percentages, with hours down the rows.
selected <- p[p$threshold %in% c(5, 10, 15, 20), ]
selected$percent <- round(100 * selected$probability, 1)
reshape(selected[c("hour_local", "threshold", "percent")],
        idvar = "hour_local", timevar = "threshold", direction = "wide")

# Dates and numeric day-of-year inputs are also supported.
tlh_wind_exceedance(as.Date("2026-06-15"), thresholds = c(5, 10, 15))
tlh_wind_exceedance(166) # non-leap reference year 2001: June 15
tlh_wind_exceedance(60, reference_year = 2000) # February 29

# Optional settings, with knots or m/s also supported.
tlh_wind_exceedance("06-15", thresholds = c(5, 10), units = "knots",
                    window_days = 2, end_year = 2025)
```

Or use the command line:

```sh
Rscript scripts/tallahassee_wind_exceedance.R 06-15 data/output/tlh_june15.csv
Rscript scripts/test_tallahassee_wind_exceedance.R
```

The command-line script saves a PNG alongside the CSV automatically. Plots
show local time (9 am–3 pm) on the horizontal axis and **P(W > threshold)**
as a percentage on the vertical axis. Each line represents a threshold
at 5 mph intervals (0, 5, 10, etc.). The CSV retains the full 1 mph table,
both probability definitions, and sample sizes. Use `thresholds = c(5, 10)`
in `plot_tlh_wind()` to select specific thresholds already in the table.
For other speed units, the default plotting interval is 5 of that unit.

**Data and method.** The station is Tallahassee airport **KTLH/TLH**,
[IEM FL_ASOS station metadata](https://mesonet.agron.iastate.edu/sites/site.php?station=TLH&network=FL_ASOS).
The script downloads routine hourly METAR sustained wind (`sknt`) from the
[Iowa Environmental Mesonet archive](https://mesonet.agron.iastate.edu/request/download.phtml),
which includes NWS/NOAA observations. The
[API documentation](https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?help=)
defines routine reports as `report_type=3`. Specials and five-minute reports
are excluded to avoid giving stormy hours more weight. These are sustained
wind samples near the hour, **not full-hour means or maximum gusts**.

- The default period is the latest **30 completed calendar years**: 1996–2025
  when run in 2026. `end_year` fixes a reproducible period. The year in a
  requested date identifies its month/day, not the historical period.
- The default seasonal window is **±2 days**, so June 15 uses June 13–17
  from every historical year. `window_days = 0` uses only the exact date.
  Each available date/hour contributes one observation with equal weight.
- UTC timestamps are converted to `America/New_York`, respecting historical
  daylight saving time. Reports are assigned to the nearest clock hour:
  08:53 represents 09:00, and 14:53 represents 15:00. If several routine
  reports map to one date/hour, use the closest valid wind report, breaking
  ties by earlier timestamp. Missing hours are not interpolated.
- For each hour and threshold `v`, the primary output `probability` is
  **P(wind <= v) = 1 - P(wind > v)**, calculated as
  `sum(wind <= v) / number_of_valid_wind_observations`. A value of 0.80 at
  10 mph means 80% of valid historical observations were at or below 10 mph.
  Equality is included. `exceedance_probability` also retains P(wind > v).
  Calm winds are retained; missing, negative, and
  nonfinite winds are excluded. Knots convert to mph using `1852/1609.344`.
- By default thresholds run from 0 to at least 30 mph in 1 mph increments,
  extending to the ceiling of the largest wind in the selected sample so
  the complete empirical cumulative probability curve reaches one. Custom thresholds
  are interpreted in the selected units. Probabilities are fractions (0–1).
- `n`, `expected_n`, `coverage`, and `years_observed` describe data support;
  `at_or_below` is the numerator for `probability`; `exceedances` counts winds
  strictly above the threshold. Most ±2-day samples have 150 possible
  observations per hour. Empty samples return `NA` probabilities, not zero.
- Windows use actual calendar days in each historical year, accounting for
  leap years and wrapping December/January. February 28 ±2 days uses
  February 26–March 2 in a non-leap year and February 26–March 1 in a leap
  year. A February 29 request uses only leap years: eight center dates or
  40 possible observations per hour with ±2 days in 1996–2025.

Raw data and a source URL/retrieval timestamp are cached in `data/cache/tlh`.
Use `refresh = TRUE` to fetch a new copy. The initial download needs internet
access; subsequent dates reuse the cache. Output CSV files include period,
units, and window settings.

These are historical climatological frequencies, not a forecast for a
particular year. Nearby days and hours can be correlated, so 150 reports
are not 150 independent weather events. A zero observed frequency does
not establish zero risk. The IEM archive has limited additional quality
control; coverage measures availability rather than measurement quality.
Changes in instruments, siting, and climate can affect a 30-year estimate.
Airport winds may differ from winds at a sheltered playing location.

## Tallahassee hourly wind quantiles

The separate quantile script uses the same station archive, 30-year period,
±2-day window, and Eastern local hours (9 am–3 pm). Run from the project root:

```r
source("scripts/tallahassee_wind_quantiles.R")
q <- tlh_wind_quantiles("10-24") # returns the table and displays the plot
q <- tlh_wind_quantiles("06-15", plot_file = "data/output/tlh_june15_quantiles.png")
plot_tlh_wind_quantiles(q) # redraw an existing table
```

Or save both a CSV and a PNG from the command line:

```sh
Rscript scripts/tallahassee_wind_quantiles.R 10-24
```

Each line shows one quantile (.90, .75, .50, .25, .10) across local hours.
The .50 line is the median. The .90 quantile is the smallest observed wind
speed for which at least 90% of valid observations are at or below that
speed. Quantiles use R's `quantile(type = 1)` (inverse empirical CDF), so
they remain observed values rather than interpolated speeds. Ties can
make the proportion at or below a quantile greater than its nominal level.
The default units are mph; `units = "knots"` or `"m/s"` is also supported.
The table includes sample sizes and coverage; missing hours return `NA`.
Use `plot = FALSE` to calculate the table without displaying a plot.
Custom `probs`, `end_year`, `window_days`, and numeric day-of-year inputs
are supported as in the probability workflow.

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
