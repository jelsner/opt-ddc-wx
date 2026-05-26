#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ncdf4)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(suncalc)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

parse_years <- function(value) {
  if (grepl(":", value, fixed = TRUE)) {
    bounds <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    if (length(bounds) != 2 || any(is.na(bounds))) {
      stop("Invalid --years value. Use a range like 2020:2025 or a list like 2020,2021.", call. = FALSE)
    }
    return(seq(bounds[[1]], bounds[[2]]))
  }

  years <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  if (any(is.na(years))) {
    stop("Invalid --years value. Use a range like 2020:2025 or a list like 2020,2021.", call. = FALSE)
  }
  years
}

parse_args <- function(args) {
  opts <- list(
    input_dir = "data/output/conus",
    years = "2020:2025",
    output_csv = "data/output/conus/conus_monthly_optimal_hours_2020_2025_average.csv",
    output_png = "data/output/conus/conus_monthly_optimal_hours_2020_2025_percentiles.png",
    percentile = 0.20,
    min_temp_f = 40,
    max_wind_mps = 5
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    name <- sub("^--", "", key)
    if (!name %in% names(opts)) {
      stop("Unknown option: ", key, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for option: ", key, call. = FALSE)
    }
    value <- args[[i + 1]]
    if (name %in% c("percentile", "min_temp_f", "max_wind_mps")) {
      value <- as.numeric(value)
    }
    opts[[name]] <- value
    i <- i + 2
  }

  opts$years <- parse_years(opts$years)
  if (is.na(opts$percentile) || opts$percentile <= 0 || opts$percentile >= 0.5) {
    stop("--percentile must be a number greater than 0 and less than 0.5.", call. = FALSE)
  }
  opts
}

get_nc_time <- function(nc) {
  time_name <- intersect(c("valid_time", "time"), names(nc$dim))[1] %||%
    intersect(c("valid_time", "time"), names(nc$var))[1]
  if (is.na(time_name)) {
    stop("Could not find a time dimension/variable named valid_time or time.", call. = FALSE)
  }

  time_raw <- ncvar_get(nc, time_name)
  units <- nc$dim[[time_name]]$units %||% nc$var[[time_name]]$units
  if (is.null(units) || grepl("^seconds since 1970-01-01", units)) {
    return(as.POSIXct(time_raw, origin = "1970-01-01", tz = "UTC"))
  }
  if (grepl("^hours since 1900-01-01", units)) {
    return(as.POSIXct(time_raw * 3600, origin = "1900-01-01", tz = "UTC"))
  }
  stop("Unsupported time units: ", units, call. = FALSE)
}

read_lat_lon <- function(nc) {
  lat_name <- intersect(c("latitude", "lat"), names(nc$dim))[1]
  lon_name <- intersect(c("longitude", "lon"), names(nc$dim))[1]
  if (is.na(lat_name) || is.na(lon_name)) {
    stop("Could not find latitude/longitude dimensions.", call. = FALSE)
  }

  lats <- ncvar_get(nc, lat_name)
  lons <- ncvar_get(nc, lon_name)
  lons <- ifelse(lons > 180, lons - 360, lons)
  list(lats = lats, lons = lons, lat_name = lat_name, lon_name = lon_name)
}

as_lat_lon_time <- function(values, var_dims, lat_name, lon_name, time_name) {
  dim_names <- vapply(var_dims, function(x) x$name, character(1))
  wanted <- c(lat_name, lon_name, time_name)
  perm <- match(wanted, dim_names)
  if (any(is.na(perm))) {
    stop("Weather variable dimensions do not include latitude, longitude, and time.", call. = FALSE)
  }
  aperm(values, perm)
}

get_weather_array <- function(nc, var_name, lat_name, lon_name, time_name) {
  values <- ncvar_get(nc, var_name, collapse_degen = FALSE)
  as_lat_lon_time(values, nc$var[[var_name]]$dim, lat_name, lon_name, time_name)
}

build_sun_cache <- function(dates, grid) {
  sun_by_date <- vector("list", length(dates))
  names(sun_by_date) <- as.character(dates)

  for (date_key in names(sun_by_date)) {
    sun_input <- grid %>% mutate(date = as.Date(date_key))
    sun_by_date[[date_key]] <- getSunlightTimes(
      data = sun_input,
      keep = c("sunrise", "sunset")
    )
  }
  sun_by_date
}

daylight_matrix <- function(datetime_utc, sun_cache, n_lat, n_lon) {
  date_utc <- as.Date(datetime_utc)
  checks <- c(as.character(date_utc), as.character(date_utc - 1))
  daylight <- rep(FALSE, n_lat * n_lon)

  for (date_key in checks) {
    sun <- sun_cache[[date_key]]
    if (!is.null(sun)) {
      daylight <- daylight | (datetime_utc >= sun$sunrise & datetime_utc <= sun$sunset)
    }
  }

  matrix(daylight, nrow = n_lat, ncol = n_lon)
}

summarise_month_file <- function(file_path, grid, opts) {
  message("Reading ", file_path)
  nc <- nc_open(file_path)
  on.exit(nc_close(nc), add = TRUE)

  coords <- read_lat_lon(nc)
  time_name <- intersect(c("valid_time", "time"), names(nc$dim))[1] %||%
    intersect(c("valid_time", "time"), names(nc$var))[1]
  timestamps <- get_nc_time(nc)

  t2m <- get_weather_array(nc, "t2m", coords$lat_name, coords$lon_name, time_name) - 273.15
  u10 <- get_weather_array(nc, "u10", coords$lat_name, coords$lon_name, time_name)
  v10 <- get_weather_array(nc, "v10", coords$lat_name, coords$lon_name, time_name)

  file_dates <- seq(min(as.Date(timestamps)) - 1, max(as.Date(timestamps)), by = "day")
  sun_cache <- build_sun_cache(file_dates, grid)
  optimal_hours <- matrix(0L, nrow = length(coords$lats), ncol = length(coords$lons))
  min_temp_c <- (opts$min_temp_f - 32) * 5 / 9

  for (t_index in seq_along(timestamps)) {
    daylight <- daylight_matrix(timestamps[[t_index]], sun_cache, length(coords$lats), length(coords$lons))
    wind_mps <- sqrt(u10[, , t_index]^2 + v10[, , t_index]^2)
    optimal_hours <- optimal_hours + (daylight & t2m[, , t_index] > min_temp_c & wind_mps < opts$max_wind_mps)
  }

  parts <- strcapture(
    "conus_era5_([0-9]{4})_([0-9]{2})\\.nc$",
    basename(file_path),
    proto = list(year = integer(), month = integer())
  )

  expand.grid(lat = coords$lats, lon = coords$lons) %>%
    mutate(
      year = parts$year,
      month = parts$month,
      optimal_hours = as.vector(optimal_hours)
    )
}

get_expected_files <- function(opts) {
  files <- file.path(
    opts$input_dir,
    sprintf("conus_era5_%d_%02d.nc", rep(opts$years, each = 12), rep(1:12, times = length(opts$years)))
  )
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Missing monthly NetCDF files:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
  files
}

average_monthly_optimal_hours <- function(opts, grid) {
  bind_rows(lapply(get_expected_files(opts), function(file) summarise_month_file(file, grid, opts))) %>%
    group_by(month, lat, lon) %>%
    summarise(
      n_years = n(),
      avg_optimal_hours = mean(optimal_hours),
      .groups = "drop"
    )
}

add_conus_percentile_class <- function(results, opts) {
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop(
      "The 'maps' R package is required to mask to the conterminous US. ",
      "Install it with install.packages('maps') and rerun this script.",
      call. = FALSE
    )
  }

  state_name <- maps::map.where("state", x = results$lon, y = results$lat)
  conus <- !is.na(state_name)

  results %>%
    mutate(conus = conus) %>%
    group_by(month) %>%
    mutate(
      lower_threshold = quantile(avg_optimal_hours[conus], opts$percentile, na.rm = TRUE, names = FALSE),
      upper_threshold = quantile(avg_optimal_hours[conus], 1 - opts$percentile, na.rm = TRUE, names = FALSE),
      percentile_class = case_when(
        conus & avg_optimal_hours <= lower_threshold ~ "Unfavorable",
        conus & avg_optimal_hours >= upper_threshold ~ "Favorable",
        TRUE ~ "Other"
      ),
      percentile_class = factor(percentile_class, levels = c("Unfavorable", "Other", "Favorable")),
      month_label = factor(month.abb[month], levels = month.abb)
    ) %>%
    ungroup()
}

filter_map_extent <- function(map_df, lon_range, lat_range, pad = 1) {
  map_df %>%
    filter(
      .data$long >= lon_range[[1]] - pad,
      .data$long <= lon_range[[2]] + pad,
      .data$lat >= lat_range[[1]] - pad,
      .data$lat <= lat_range[[2]] + pad
    )
}

get_border_data <- function(lon_range, lat_range) {
  list(
    states = filter_map_extent(ggplot2::map_data("state"), lon_range, lat_range),
    countries = filter_map_extent(
      ggplot2::map_data("world", region = c("USA", "Canada", "Mexico")),
      lon_range,
      lat_range
    )
  )
}

plot_monthly_percentiles <- function(results) {
  lon_range <- range(results$lon)
  lat_range <- range(results$lat)
  borders <- get_border_data(lon_range, lat_range)
  caption <- paste(
    "Based on the number of daylight hours (2020-2025) warmer than 4 C (40 F)",
    "and with winds lighter than 5 m/s (11 mph).",
    "ERA5 data from climate.copernicus.eu/climate-reanalysis"
  )

  ggplot(results, aes(x = lon, y = lat, fill = percentile_class)) +
    geom_raster() +
    geom_polygon(
      data = borders$states,
      aes(x = long, y = lat, group = group),
      inherit.aes = FALSE,
      fill = NA,
      color = "grey45",
      linewidth = 0.15
    ) +
    geom_polygon(
      data = borders$countries,
      aes(x = long, y = lat, group = group),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.25
    ) +
    facet_wrap(~month_label, ncol = 4) +
    coord_quickmap(xlim = lon_range, ylim = lat_range, expand = FALSE) +
    scale_fill_manual(
      name = NULL,
      values = c("Unfavorable" = "#8b5a2b", "Other" = "grey82", "Favorable" = "#2f855a"),
      breaks = c("Unfavorable", "Favorable"),
      labels = c("Least favorable (Bottom 20%)", "Most favorable (Top 20%)"),
      drop = FALSE
    ) +
    labs(
      title = "Favorable and unfavorable places to play DDC in the continental U.S.",
      x = NULL,
      y = NULL,
      caption = caption
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 16),
      plot.caption = element_text(hjust = 0, size = 9),
      strip.text = element_text(face = "bold")
    )
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  first_nc <- nc_open(file.path(opts$input_dir, sprintf("conus_era5_%d_01.nc", opts$years[[1]])))
  coords <- read_lat_lon(first_nc)
  nc_close(first_nc)
  grid <- expand.grid(lat = coords$lats, lon = coords$lons)

  results <- average_monthly_optimal_hours(opts, grid) %>%
    add_conus_percentile_class(opts)

  dir.create(dirname(opts$output_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(results, opts$output_csv, row.names = FALSE)

  p <- plot_monthly_percentiles(results)
  ggsave(opts$output_png, p, width = 14, height = 10, dpi = 300, bg = "white")

  message("Wrote CSV: ", opts$output_csv)
  message("Wrote plot: ", opts$output_png)
}

main()
