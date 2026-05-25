#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
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
    output_csv = "data/output/conus/conus_optimal_hours_2020_2025_average.csv",
    output_png = "data/output/conus/conus_optimal_hours_2020_2025_percentiles.png",
    percentile = 0.10,
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
    if (name == "percentile") {
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

read_year_file <- function(input_dir, year) {
  file_path <- file.path(input_dir, paste0("conus_optimal_hours_", year, ".csv"))
  if (!file.exists(file_path)) {
    stop("Missing annual summary CSV: ", file_path, call. = FALSE)
  }

  df <- read.csv(file_path)
  required <- c("lat", "lon", "optimal_hours")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }

  df %>%
    select(lat, lon, optimal_hours) %>%
    mutate(year = year)
}

average_optimal_hours <- function(opts) {
  bind_rows(lapply(opts$years, function(year) read_year_file(opts$input_dir, year))) %>%
    group_by(lat, lon) %>%
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
  thresholds <- quantile(
    results$avg_optimal_hours[conus],
    probs = c(opts$percentile, 1 - opts$percentile),
    na.rm = TRUE,
    names = FALSE
  )

  results %>%
    mutate(
      conus = conus,
      percentile_class = case_when(
        conus & avg_optimal_hours <= thresholds[[1]] ~ "Lower 10%",
        conus & avg_optimal_hours >= thresholds[[2]] ~ "Upper 10%",
        TRUE ~ "Other"
      ),
      percentile_class = factor(percentile_class, levels = c("Lower 10%", "Other", "Upper 10%")),
      lower_threshold = thresholds[[1]],
      upper_threshold = thresholds[[2]]
    )
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

plot_percentiles <- function(results, opts) {
  lon_range <- range(results$lon)
  lat_range <- range(results$lat)
  borders <- get_border_data(lon_range, lat_range)
  year_label <- paste0(min(opts$years), "-", max(opts$years))
  title <- opts$title %||% paste0(
    "Average Optimal Hours: Upper and Lower ",
    opts$percentile * 100,
    "% of CONUS Grid Cells (",
    year_label,
    ")"
  )

  ggplot(results, aes(x = lon, y = lat, fill = percentile_class)) +
    geom_raster() +
    geom_polygon(
      data = borders$states,
      aes(x = long, y = lat, group = group),
      inherit.aes = FALSE,
      fill = NA,
      color = "grey35",
      linewidth = 0.25
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
    scale_fill_manual(
      name = NULL,
      values = c("Lower 10%" = "#2b6cb0", "Other" = "grey82", "Upper 10%" = "#c53030"),
      drop = FALSE
    ) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  results <- average_optimal_hours(opts) %>%
    add_conus_percentile_class(opts)

  dir.create(dirname(opts$output_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(results, opts$output_csv, row.names = FALSE)

  p <- plot_percentiles(results, opts)
  ggsave(opts$output_png, p, width = 11, height = 7, dpi = 300, bg = "white")

  message("Wrote CSV: ", opts$output_csv)
  message("Wrote plot: ", opts$output_png)
}

main()
