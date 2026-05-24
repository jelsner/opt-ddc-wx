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

parse_args <- function(args) {
  opts <- list(
    input_dir = "data/output/conus",
    output_csv = "data/output/conus/conus_optimal_hours.csv",
    output_png = "data/output/conus/conus_optimal_hours.png",
    pattern = "^conus_era5_\\d{4}_\\d{2}\\.nc$",
    min_temp_f = -60,
    max_wind_mps = 5,
    title = NULL
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
    if (name %in% c("min_temp_f", "max_wind_mps")) {
      value <- as.numeric(value)
    }
    opts[[name]] <- value
    i <- i + 2
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
  if (is.null(units)) {
    return(as.POSIXct(time_raw, origin = "1970-01-01", tz = "UTC"))
  }

  if (grepl("^seconds since 1970-01-01", units)) {
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

summarise_file <- function(file_path, grid, opts, totals) {
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

  for (t_index in seq_along(timestamps)) {
    daylight <- daylight_matrix(
      timestamps[[t_index]],
      sun_cache,
      length(coords$lats),
      length(coords$lons)
    )
    wind_mps <- sqrt(u10[, , t_index]^2 + v10[, , t_index]^2)
    optimal <- daylight & t2m[, , t_index] > (opts$min_temp_f - 32) * 5 / 9 & wind_mps < opts$max_wind_mps

    totals$daylight <- totals$daylight + daylight
    totals$optimal <- totals$optimal + optimal
  }

  totals
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
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop(
      "The 'maps' R package is required to draw state, county, and country borders. ",
      "Install it with install.packages('maps') and rerun this script.",
      call. = FALSE
    )
  }

  list(
    counties = filter_map_extent(ggplot2::map_data("county"), lon_range, lat_range),
    states = filter_map_extent(ggplot2::map_data("state"), lon_range, lat_range),
    countries = filter_map_extent(
      ggplot2::map_data("world", region = c("USA", "Canada", "Mexico")),
      lon_range,
      lat_range
    )
  )
}

plot_results <- function(results, opts) {
  title <- opts$title %||% paste0(
    "Optimal Playable Daylight Hours: temp > ",
    opts$min_temp_f,
    " F and wind < ",
    opts$max_wind_mps,
    " m/s"
  )
  lon_range <- range(results$lon)
  lat_range <- range(results$lat)
  if (diff(lon_range) == 0) lon_range <- lon_range + c(-0.25, 0.25)
  if (diff(lat_range) == 0) lat_range <- lat_range + c(-0.25, 0.25)
  borders <- get_border_data(lon_range, lat_range)

  ggplot(results, aes(x = lon, y = lat, fill = optimal_hours)) +
    geom_raster() +
    geom_polygon(
      data = borders$counties,
      aes(x = long, y = lat, group = group),
      inherit.aes = FALSE,
      fill = NA,
      color = "grey35",
      linewidth = 0.08,
      alpha = 0.35
    ) +
    geom_polygon(
      data = borders$states,
      aes(x = long, y = lat, group = group),
      inherit.aes = FALSE,
      fill = NA,
      color = "grey15",
      linewidth = 0.25,
      alpha = 0.75
    ) +
    geom_polygon(
      data = borders$countries,
      aes(x = long, y = lat, group = group),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.45
    ) +
    coord_quickmap(xlim = lon_range, ylim = lat_range, expand = FALSE) +
    scale_fill_viridis_c(name = "Hours", option = "C", na.value = "grey90") +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  files <- list.files(opts$input_dir, pattern = opts$pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No NetCDF files matched ", opts$pattern, " in ", opts$input_dir, call. = FALSE)
  }

  first_nc <- nc_open(files[[1]])
  coords <- read_lat_lon(first_nc)
  nc_close(first_nc)
  grid <- expand.grid(lat = coords$lats, lon = coords$lons)
  totals <- list(
    daylight = matrix(0L, nrow = length(coords$lats), ncol = length(coords$lons)),
    optimal = matrix(0L, nrow = length(coords$lats), ncol = length(coords$lons))
  )

  for (file in sort(files)) {
    totals <- summarise_file(file, grid, opts, totals)
  }

  results <- expand.grid(lat = coords$lats, lon = coords$lons) %>%
    mutate(
      daylight_hours = as.vector(totals$daylight),
      optimal_hours = as.vector(totals$optimal),
      optimal_fraction = if_else(daylight_hours > 0, optimal_hours / daylight_hours, NA_real_)
    )

  dir.create(dirname(opts$output_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(results, opts$output_csv, row.names = FALSE)

  p <- plot_results(results, opts)
  ggsave(opts$output_png, p, width = 11, height = 7, dpi = 300)

  message("Wrote CSV: ", opts$output_csv)
  message("Wrote plot: ", opts$output_png)
}

main()
