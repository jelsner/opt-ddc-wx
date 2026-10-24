# Base R only. Reuses the station cache and preparation from the sibling script.
.tlh_quantile_source <- if (sys.nframe() > 0L && !is.null(sys.frame(1)$ofile)) {
  sys.frame(1)$ofile
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
}
source(file.path(dirname(normalizePath(.tlh_quantile_source)),
                 "tallahassee_wind_exceedance.R"), local = TRUE)
rm(.tlh_quantile_source)

tlh_wind_quantiles <- function(day, probs = c(.9, .75, .5, .25, .1),
                               units = c("mph", "knots", "m/s"),
                               window_days = 2L, end_year = NULL,
                               reference_year = 2001L,
                               cache_dir = "data/cache/tlh", refresh = FALSE,
                               observations = NULL, plot = TRUE,
                               plot_file = NULL) {
  units <- match.arg(units)
  md <- tlh_day(day, reference_year)
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs) | probs < 0 | probs > 1))
    stop("probs must contain quantile probabilities from 0 to 1.")
  if (length(window_days) != 1L || !is.finite(window_days) ||
      window_days < 0 || window_days > 182 || window_days != floor(window_days))
    stop("window_days must be an integer from 0 to 182.")
  if (is.null(end_year)) end_year <- max(tlh_years())
  if (length(end_year) != 1L || !is.finite(end_year) || end_year != floor(end_year) ||
      end_year < 1937 || end_year > max(tlh_years()))
    stop("end_year must be a completed year >= 1937.")
  years <- seq.int(end_year - 29L, end_year)
  if (is.null(observations)) observations <- tlh_download(years, cache_dir, refresh)
  x <- tlh_prepare(observations, years)
  x <- x[tlh_select_dates(x$date, md, window_days), ]
  dates <- seq(as.Date(sprintf("%d-01-01", min(years))),
               as.Date(sprintf("%d-12-31", max(years))), by = "day")
  expected <- sum(tlh_select_dates(dates, md, window_days))
  factor <- switch(units, mph = 1852 / 1609.344, knots = 1, "m/s" = 1852 / 3600)
  result <- do.call(rbind, lapply(9:15, function(h) {
    a <- x[x$hour == h & !is.na(x$wind_knots), ]
    n <- nrow(a)
    # Type 1 is the inverse empirical CDF: quantiles remain observed speeds.
    wind <- if (n) as.numeric(stats::quantile(a$wind_knots * factor,
                        probs = probs, type = 1, names = FALSE)) else rep(NA_real_, length(probs))
    data.frame(month_day = md, hour_local = h, time_zone = "America/New_York",
               quantile = probs, wind_speed = wind, units = units,
               n = n, expected_n = expected, coverage = n / expected,
               years_observed = length(unique(a$year)), start_year = min(years),
               end_year = max(years), window_days = window_days, quantile_type = 1L)
  }))
  if (any(result$coverage < .9, na.rm = TRUE))
    warning("Some hours have less than 90% coverage; inspect n and years_observed.")
  if (any(result$n == 0L)) warning("Some hours have no valid observations.")
  rownames(result) <- NULL
  if (plot || !is.null(plot_file)) plot_tlh_wind_quantiles(result, file = plot_file)
  result
}

plot_tlh_wind_quantiles <- function(quantiles, file = NULL) {
  q <- quantiles
  required <- c("month_day", "hour_local", "quantile", "wind_speed", "units",
                "start_year", "end_year", "window_days")
  if (!all(required %in% names(q)) || !nrow(q)) stop("Supply an hourly wind quantile table.")
  for (name in c("month_day", "units", "start_year", "end_year", "window_days")) {
    if (length(unique(q[[name]])) != 1L) stop("Plot one date and set of settings at a time.")
  }
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = 1600, height = 1000, res = 160)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old <- graphics::par(mar = c(5, 5, 5, 9), xpd = FALSE)
  on.exit(graphics::par(old), add = TRUE, after = FALSE)
  probs <- sort(unique(q$quantile), decreasing = TRUE)
  colors <- grDevices::hcl.colors(length(probs), "Dark 3")
  hours <- sort(unique(q$hour_local))
  upper <- max(c(5, q$wind_speed), na.rm = TRUE)
  graphics::plot(c(9, 15), c(0, ceiling(upper / 5) * 5), type = "n", xaxt = "n",
                 xlab = "Local time (Eastern)",
                 ylab = paste0("Sustained wind speed (", q$units[1], ")"),
                 main = paste0("Tallahassee hourly wind quantiles | ", q$month_day[1]),
                 las = 1)
  graphics::axis(1, at = hours, labels = sprintf("%d %s", (hours + 11L) %% 12L + 1L,
                        ifelse(hours < 12L, "am", "pm")))
  graphics::abline(h = graphics::axTicks(2), col = "#E4E8EC")
  graphics::abline(v = hours, col = "#F0F2F4")
  for (i in seq_along(probs)) {
    a <- q[q$quantile == probs[i], ]
    a <- a[order(a$hour_local), ]
    graphics::lines(a$hour_local, a$wind_speed, col = colors[i],
                    lwd = if (probs[i] == .5) 3 else 2, type = "o", pch = 14 + i, cex = .8)
  }
  graphics::mtext(sprintf("KTLH airport | %d–%d | ±%d calendar days | Eastern local time",
                         q$start_year[1], q$end_year[1], q$window_days[1]),
                  side = 3, line = .6, cex = .85)
  bounds <- graphics::par("usr")
  graphics::legend(x = bounds[2] + .025 * diff(bounds[1:2]), y = bounds[4],
                    xjust = 0, yjust = 1, xpd = TRUE, bty = "n", cex = .85,
                    legend = paste0(100 * probs, "%"), title = "Quantile",
                    col = colors, lwd = ifelse(probs == .5, 3, 2),
                    pch = 14 + seq_along(probs))
  invisible(q)
}

# Rscript scripts/tallahassee_wind_quantiles.R MM-DD [output.csv]
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Usage: Rscript scripts/tallahassee_wind_quantiles.R MM-DD [output.csv]")
  output <- if (length(args) >= 2L) args[2] else paste0("data/output/tlh_quantiles_", args[1], ".csv")
  image_path <- paste0(tools::file_path_sans_ext(output), ".png")
  answer <- tlh_wind_quantiles(args[1], plot_file = image_path)
  write.csv(answer, output, row.names = FALSE)
  print(answer, row.names = FALSE)
  message("Saved ", output, " and ", image_path)
}
