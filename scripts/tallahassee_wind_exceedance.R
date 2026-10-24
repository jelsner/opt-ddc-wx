# Base R only. Source this file, then call tlh_wind_exceedance("06-15").
# Routine METAR sustained wind, not a full-hour mean or gust maximum.
# Main output probability = P(wind <= threshold) = 1 - P(wind > threshold).

tlh_years <- function(as_of = Sys.Date()) {
  last <- as.integer(format(as.Date(as_of), "%Y")) - 1L
  seq.int(last - 29L, last)
}

tlh_url <- function(years) {
  paste0("https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?",
         "station=TLH&network=FL_ASOS&data=sknt&report_type=3",
         "&tz=UTC&format=onlycomma&latlon=no&elev=no&missing=M",
         "&year1=", min(years), "&month1=1&day1=1",
         "&year2=", max(years) + 1L, "&month2=1&day2=1")
}

tlh_read_raw <- function(path) {
  x <- read.csv(path, na.strings = c("M", "null", ""),
                stringsAsFactors = FALSE, comment.char = "#")
  if (!all(c("station", "valid", "sknt") %in% names(x)) || !nrow(x))
    stop("Invalid or empty IEM response: ", path)
  if (any(is.na(x$station) | x$station != "TLH"))
    stop("Unexpected station in IEM response.")
  x
}

tlh_download <- function(years = tlh_years(), cache_dir = "data/cache/tlh",
                         refresh = FALSE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(cache_dir, sprintf("TLH_routine_UTC_%d_%d.csv",
                                      min(years), max(years)))
  if (!file.exists(path) || refresh) {
    tmp <- tempfile(tmpdir = cache_dir, fileext = ".csv")
    on.exit(unlink(tmp), add = TRUE)
    old <- options(timeout = max(600, getOption("timeout")))
    on.exit(options(old), add = TRUE)
    url <- tlh_url(years)
    message("Downloading TLH routine reports for ", min(years), "–", max(years))
    status <- download.file(url, tmp, mode = "wb", quiet = TRUE)
    if (status != 0L) stop("IEM download failed; retry later.")
    tlh_read_raw(tmp) # Validate before replacing the cache.
    if (!file.copy(tmp, path, overwrite = TRUE)) stop("Cannot save cache.")
    writeLines(c(paste("URL:", url), paste("Retrieved UTC:",
               format(Sys.time(), tz = "UTC", usetz = TRUE))),
               paste0(path, ".source.txt"))
  }
  tlh_read_raw(path)
}

tlh_day <- function(day, reference_year = 2001L) {
  if (inherits(day, "Date")) {
    md <- format(day, "%m-%d")
  } else if (is.numeric(day)) {
    max_day <- as.integer(as.Date(sprintf("%d-01-01", reference_year + 1L)) -
                            as.Date(sprintf("%d-01-01", reference_year)))
    if (length(day) != 1L || !is.finite(day) || day != floor(day) ||
        day < 1 || day > max_day) stop("Invalid numeric day of year.")
    md <- format(as.Date(sprintf("%d-01-01", reference_year)) + day - 1, "%m-%d")
  } else {
    md <- as.character(day)
    if (length(md) == 1L && grepl("^\\d{4}-\\d{2}-\\d{2}$", md)) {
      d <- as.Date(md, format = "%Y-%m-%d")
      if (is.na(d) || format(d, "%Y-%m-%d") != md) stop("Invalid date.")
      md <- format(d, "%m-%d")
    }
  }
  if (length(md) != 1L || is.na(md) || !grepl("^\\d{2}-\\d{2}$", md))
    stop("Use MM-DD, YYYY-MM-DD, a Date, or a numeric day of year.")
  d <- as.Date(paste0("2000-", md), format = "%Y-%m-%d")
  if (is.na(d) || format(d, "%m-%d") != md) stop("Invalid calendar day.")
  md
}

tlh_prepare <- function(raw, years) {
  time <- as.POSIXct(raw$valid, format = "%Y-%m-%d %H:%M", tz = "UTC")
  if (anyNA(time)) stop("Unparseable IEM timestamps.")
  # Round to the nearest clock hour; e.g., 08:53 represents 09:00.
  clock <- as.POSIXct(floor((as.numeric(time) + 1800) / 3600) * 3600,
                     origin = "1970-01-01", tz = "UTC")
  local <- as.POSIXlt(clock, tz = "America/New_York")
  wind <- suppressWarnings(as.numeric(raw$sknt))
  wind[!is.finite(wind) | wind < 0] <- NA_real_
  x <- data.frame(date = as.Date(format(clock, "%Y-%m-%d", tz = "America/New_York")),
                  year = local$year + 1900L, hour = local$hour,
                  offset = abs(as.numeric(time) - as.numeric(clock)),
                  wind_knots = wind, time_utc = time)
  x <- x[x$year %in% years & x$hour %in% 9:15, ]
  # Prefer a valid report, then the closest report; never count a date/hour twice.
  x <- x[order(x$date, x$hour, is.na(x$wind_knots), x$offset, x$time_utc), ]
  x <- x[!duplicated(x[c("date", "hour")]), ]
  x$month_day <- format(x$date, "%m-%d")
  x
}

tlh_select_dates <- function(dates, month_day, window_days) {
  year <- as.integer(format(dates, "%Y"))
  target <- as.Date(paste0(year, "-", month_day), format = "%Y-%m-%d")
  year_length <- as.integer(as.Date(paste0(year + 1L, "-01-01")) -
                              as.Date(paste0(year, "-01-01")))
  delta <- abs(as.integer(dates - target))
  # Circular distance wraps New Year. Feb 29 has no target in non-leap years.
  !is.na(target) & pmin(delta, year_length - delta) <= window_days
}

tlh_wind_exceedance <- function(day, thresholds = NULL,
                                units = c("mph", "knots", "m/s"),
                                window_days = 2L, end_year = NULL,
                                reference_year = 2001L,
                                cache_dir = "data/cache/tlh", refresh = FALSE,
                                observations = NULL, plot = FALSE,
                                plot_file = NULL) {
  units <- match.arg(units)
  md <- tlh_day(day, reference_year)
  if (!is.null(thresholds) && (!is.numeric(thresholds) || !length(thresholds) ||
      any(!is.finite(thresholds) | thresholds < 0))) stop("Invalid thresholds.")
  if (length(window_days) != 1L || !is.finite(window_days) ||
      window_days < 0 || window_days > 182 || window_days != floor(window_days))
    stop("window_days must be an integer from 0 to 182.")
  if (is.null(end_year)) end_year <- max(tlh_years())
  if (length(end_year) != 1L || !is.finite(end_year) ||
      end_year != floor(end_year) || end_year < 1937 ||
      end_year > max(tlh_years())) stop("end_year must be a completed year >= 1937.")
  years <- seq.int(end_year - 29L, end_year)
  if (is.null(observations)) observations <- tlh_download(years, cache_dir, refresh)
  x <- tlh_prepare(observations, years)
  select_dates <- function(dates) {
    tlh_select_dates(dates, md, window_days)
  }
  x <- x[select_dates(x$date), ]
  dates <- seq(as.Date(sprintf("%d-01-01", min(years))),
               as.Date(sprintf("%d-12-31", max(years))), by = "day")
  expected <- sum(select_dates(dates))
  factor <- switch(units, mph = 1852 / 1609.344, knots = 1, "m/s" = 1852 / 3600)
  if (is.null(thresholds)) {
    speeds <- x$wind_knots[!is.na(x$wind_knots)] * factor
    thresholds <- seq.int(0L, ceiling(max(c(30, speeds))))
  }
  result <- do.call(rbind, lapply(9:15, function(h) {
    a <- x[x$hour == h & !is.na(x$wind_knots), ]
    n <- nrow(a)
    k <- vapply(thresholds, function(v) sum(a$wind_knots * factor > v), integer(1))
    data.frame(month_day = md, hour_local = h,
               time_zone = "America/New_York", threshold = thresholds, units = units,
               probability = if (n) 1 - k / n else NA_real_,
               exceedance_probability = if (n) k / n else NA_real_,
               at_or_below = n - k, exceedances = k,
               n = n, expected_n = expected,
               coverage = if (expected) n / expected else NA_real_,
               years_observed = length(unique(a$year)),
               start_year = min(years), end_year = max(years),
               window_days = window_days)
  }))
  if (any(result$coverage < 0.9, na.rm = TRUE))
    warning("Some hours have less than 90% coverage; inspect n and years_observed.")
  if (any(result$n == 0)) warning("Some hours have no valid observations.")
  rownames(result) <- NULL
  if (plot || !is.null(plot_file)) plot_tlh_wind(result, file = plot_file)
  result
}

# Hourly probability curves at thresholds spaced every 5 speed units.
plot_tlh_wind <- function(probabilities, file = NULL,
                          type = c("exceedance", "at_or_below"),
                          thresholds = NULL) {
  type <- match.arg(type)
  p <- probabilities
  required <- c("month_day", "hour_local", "units", "threshold", "probability",
                "exceedance_probability", "start_year", "end_year", "window_days", "n")
  if (!all(required %in% names(p)) || !nrow(p)) stop("Supply a wind probability table.")
  for (name in c("month_day", "units", "start_year", "end_year", "window_days")) {
    if (length(unique(p[[name]])) != 1L) stop("Plot one date and set of settings at a time.")
  }
  if (is.null(thresholds))
    thresholds <- sort(unique(p$threshold[abs(p$threshold / 5 - round(p$threshold / 5)) < 1e-8]))
  if (!length(thresholds) || anyNA(thresholds) || !all(thresholds %in% p$threshold))
    stop("Plot thresholds must be present in the probability table.")
  thresholds <- sort(unique(thresholds))
  p <- p[p$threshold %in% thresholds, ]
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = 1600, height = 1000, res = 160)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old <- graphics::par(mar = c(5, 5, 5, 9), xpd = FALSE)
  on.exit(graphics::par(old), add = TRUE, after = FALSE)
  hours <- sort(unique(p$hour_local))
  colors <- grDevices::hcl.colors(length(thresholds), "Dark 3")
  column <- if (type == "at_or_below") "probability" else "exceedance_probability"
  description <- if (type == "at_or_below") "P(W <= threshold)" else "P(W > threshold)"
  graphics::plot(c(9, 15), c(0, 100), type = "n", yaxt = "n", xaxt = "n",
                 xlab = "Local time (Eastern)",
                 ylab = paste0(description, " (%)"),
                 main = paste0("Tallahassee hourly wind probabilities | ", p$month_day[1]),
                 sub = "")
  hour_labels <- sprintf("%d %s", (hours + 11L) %% 12L + 1L,
                         ifelse(hours < 12L, "am", "pm"))
  graphics::axis(1, at = hours, labels = hour_labels)
  graphics::axis(2, at = seq(0, 100, 20), las = 1)
  graphics::abline(h = seq(0, 100, 20), col = "#E4E8EC")
  graphics::abline(v = hours, col = "#F0F2F4")
  labels <- paste(thresholds, p$units[1])
  for (i in seq_along(thresholds)) {
    a <- p[p$threshold == thresholds[i], ]
    a <- a[order(a$hour_local), ]
    graphics::lines(a$hour_local, 100 * a[[column]], col = colors[i], lwd = 2,
                    lty = 1 + (i - 1L) %/% 4L, type = "o", pch = 15 + (i - 1L) %% 5L,
                    cex = 0.8)
  }
  graphics::mtext(sprintf("KTLH airport | %d–%d | ±%d calendar days | Eastern local time",
                         p$start_year[1], p$end_year[1], p$window_days[1]),
                  side = 3, line = 0.6, cex = 0.85)
  bounds <- graphics::par("usr")
  graphics::legend(x = bounds[2] + 0.025 * (bounds[2] - bounds[1]),
                    y = bounds[4], xjust = 0, yjust = 1, xpd = TRUE,
                    legend = labels, col = colors, lwd = 2,
                    lty = 1 + (seq_along(thresholds) - 1L) %/% 4L,
                    pch = 15 + (seq_along(thresholds) - 1L) %% 5L,
                    title = "Threshold", bty = "n", cex = 0.8)
  invisible(p)
}

# Rscript scripts/tallahassee_wind_exceedance.R MM-DD [output.csv]
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Usage: Rscript scripts/tallahassee_wind_exceedance.R MM-DD [output.csv]")
  answer <- tlh_wind_exceedance(args[1])
  output <- if (length(args) >= 2L) args[2] else paste0("data/output/tlh_wind_", args[1], ".csv")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  write.csv(answer, output, row.names = FALSE)
  image_path <- paste0(tools::file_path_sans_ext(output), ".png")
  plot_tlh_wind(answer, file = image_path)
  print(answer, row.names = FALSE)
  message("Saved ", output, " and ", image_path)
}
