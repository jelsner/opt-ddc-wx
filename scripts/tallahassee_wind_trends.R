# Quantile regression trends at 09:00–15:00 Eastern, using 30 completed years.
# Run from the project root: Rscript scripts/tallahassee_wind_trends.R [MM-DD] [B]
source("scripts/tallahassee_wind_exceedance.R")
if (dir.exists("data/r-library")) .libPaths(c("data/r-library", .libPaths()))

tlh_wind_trends <- function(day = "10-24", window_days = 2L,
                            end_year = max(tlh_years()), B = 999L, seed = 20261025L,
                            output_dir = "data/output/tlh_trends", cores = 4L) {
  if (!requireNamespace("quantreg", quietly = TRUE))
    stop("Install the quantreg package before running this script.")
  if (B < 199L || B != as.integer(B)) stop("Use at least 199 bootstrap replicates.")
  md <- tlh_day(day)
  years <- seq.int(end_year - 29L, end_year)
  raw <- tlh_download(years)
  x <- tlh_prepare(raw, years)
  x$wind_mph <- x$wind_knots * 1852 / 1609.344
  # Change in mph per decade; year is centered at the period midpoint.
  x$decade <- (x$year - mean(range(years))) / 10
  doy <- as.integer(format(x$date, "%j")) - 1
  year_length <- as.integer(as.Date(paste0(x$year + 1L, "-01-01")) -
                              as.Date(paste0(x$year, "-01-01")))
  phase <- 2 * pi * doy / year_length
  x$sin1 <- sin(phase)
  x$cos1 <- cos(phase)
  x$sin2 <- sin(2 * phase)
  x$cos2 <- cos(2 * phase)
  x <- x[!is.na(x$wind_mph), ]
  scopes <- list(all_dates = x,
                 date_window = x[tlh_select_dates(x$date, md, window_days), ])
  jobs <- expand.grid(tau = c(.5, .75, .9), h = 9:15, scope = names(scopes),
                      stringsAsFactors = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  worker <- function(index) {
    set.seed(seed + index)
    scope <- jobs$scope[index]
    h <- jobs$h[index]
    tau <- jobs$tau[index]
    warnings <- character()
    data <- scopes[[scope]]
      a <- data[data$hour == h, ]
      groups <- length(unique(a$year))
      if (groups < 20L) stop("Too few years for cluster inference: ", scope, ", hour ", h)
      formula <- if (scope == "all_dates")
        wind_mph ~ decade + sin1 + cos1 + sin2 + cos2 else wind_mph ~ decade
        key <- paste(scope, h, tau, sep = "_")
        checkpoint <- file.path(output_dir, paste0("fit_", key, ".rds"))
        if (file.exists(checkpoint)) {
          cached <- readRDS(checkpoint)
          if (identical(cached$settings, list(md, window_days, end_year, B, seed + index)))
            return(cached)
        }
        message("Fitting ", key, " (", nrow(a), " observations; ", B, " bootstrap draws)")
        fit_and_summary <- withCallingHandlers({
          fit <- quantreg::rq(formula, data = a, tau = tau, method = "fn")
          inference <- summary(fit, se = "boot", cluster = a$year, R = B, covariance = TRUE)
          list(fit = fit, inference = inference)
        }, warning = function(w) {
          warnings <<- c(warnings, paste(key, conditionMessage(w), sep = ": "))
          invokeRestart("muffleWarning")
        })
        fit <- fit_and_summary$fit
        summary <- fit_and_summary$inference
        slope <- unname(coef(fit)["decade"])
        se <- unname(summary$coefficients["decade", 2])
        # Bootstrap SE with a t reference based on independent year clusters.
        critical <- qt(.975, df = groups - 1L)
        p <- 2 * pt(-abs(slope / se), df = groups - 1L)
        result <- data.frame(scope = scope, month_day = if (scope == "date_window") md else NA,
          window_days = if (scope == "date_window") window_days else NA,
          hour_local = h, quantile = tau, slope_mph_decade = slope,
          se_mph_decade = se, ci_low = slope - critical * se,
          ci_high = slope + critical * se, p_value = p,
          n = nrow(a), year_clusters = groups, start_year = min(years), end_year = max(years),
          fitted_change_mph = slope * diff(range(years)) / 10,
          bootstrap_replicates = B, seed = seed + index)
        fit_details <- list(coefficients = coef(fit), bootstrap_coefficients = summary$B,
                            formula = deparse(formula), warnings = warnings[startsWith(warnings, key)])
        value <- list(result = result, fit = fit_details, warnings = warnings,
                       settings = list(md, window_days, end_year, B, seed + index))
        saveRDS(value, checkpoint)
        value
  }
  cores <- if (.Platform$OS.type == "windows") 1L else max(1L, min(as.integer(cores), 4L))
  completed <- parallel::mclapply(seq_len(nrow(jobs)), worker, mc.cores = cores,
                                  mc.preschedule = FALSE, mc.set.seed = FALSE)
  if (any(vapply(completed, inherits, logical(1), "try-error")))
    stop("A model failed; completed models were checkpointed. Inspect warnings and retry.")
  out <- do.call(rbind, lapply(completed, `[[`, "result"))
  fits <- lapply(completed, `[[`, "fit")
  names(fits) <- paste(jobs$scope, jobs$h, jobs$tau, sep = "_")
  warnings <- unlist(lapply(completed, `[[`, "warnings"))
  rownames(out) <- NULL
  out$p_bh_within_scope <- ave(out$p_value, out$scope, FUN = function(p) p.adjust(p, "BH"))
  out$p_bh_all <- p.adjust(out$p_value, "BH")
  out$significant_05 <- out$p_bh_all < .05
  out$direction <- ifelse(abs(out$slope_mph_decade) < 1e-7, "flat",
                           ifelse(out$slope_mph_decade > 0, "upward", "downward"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(out, file.path(output_dir, "quantile_trends.csv"), row.names = FALSE)
  saveRDS(fits, file.path(output_dir, "quantile_trend_fits.rds"))
  coverage <- aggregate(wind_mph ~ year + hour, x, length)
  names(coverage)[3] <- "valid_n"
  write.csv(coverage, file.path(output_dir, "annual_hourly_coverage.csv"), row.names = FALSE)
  writeLines(c(paste("quantreg", utils::packageVersion("quantreg")),
               "Cluster wild-gradient bootstrap SE; year clusters; t(G-1) p-values and 95% CIs.",
               "Benjamini-Hochberg adjustment across all 42 tests is the primary significance criterion.",
               "All-date models include annual and semiannual sine/cosine terms.",
               "Date-window models have an intercept and linear year trend only.",
               "Pointwise CIs are unadjusted for multiple testing.",
               "Models assume a constant linear trend; station changes can confound climate trends.",
               "Warnings:", unique(warnings), capture.output(sessionInfo())),
               file.path(output_dir, "analysis_notes.txt"))
  for (scope in names(scopes)) plot_tlh_wind_trends(out[out$scope == scope, ],
                                    file.path(output_dir, paste0(scope, "_slopes.png")))
  tlh_trend_report(out, file.path(output_dir, "report.md"))
  out
}

tlh_trend_report <- function(results, file) {
  lines <- c("# Tallahassee hourly wind-speed trends", "",
    sprintf("Period: %d–%d. Station: TLH/KTLH airport. Hours: 9 am–3 pm Eastern.",
            results$start_year[1], results$end_year[1]), "",
    sprintf("**%d of %d tests are significant after Benjamini–Hochberg adjustment at 5%%.**",
            sum(results$significant_05, na.rm = TRUE), nrow(results)), "",
    "The two scopes are all dates (adjusted for seasonality) and the date window used in the earlier plots.",
    "Slopes are mph per decade. A positive slope indicates increasing wind speed at that quantile.", "",
    "| Scope | Quantile | Slope range across hours | Raw p < .05 | Adjusted p < .05 |",
    "|---|---:|---:|---:|---:|")
  for (scope in unique(results$scope)) for (tau in c(.5, .75, .9)) {
    a <- results[results$scope == scope & results$quantile == tau, ]
    label <- if (scope == "all_dates") "All dates" else paste0(a$month_day[1], " ±", a$window_days[1], " days")
    lines <- c(lines, sprintf("| %s | %.0f%% | %.3f to %.3f | %d / 7 | %d / 7 |",
      label, 100 * tau, min(a$slope_mph_decade), max(a$slope_mph_decade),
      sum(a$p_value < .05, na.rm = TRUE), sum(a$significant_05, na.rm = TRUE)))
  }
  for (scope in unique(results$scope)) {
    a <- results[results$scope == scope, ]
    lines <- c(lines, "", paste0("## ", if (scope == "all_dates") "All dates" else
                    paste0(a$month_day[1], " ±", a$window_days[1], " days")), "",
      "Each entry is slope [pointwise 95% confidence interval], adjusted p-value.", "",
      "| Local hour | 50th percentile | 75th percentile | 90th percentile |",
      "|---|---|---|---|")
    for (h in 9:15) {
      entries <- vapply(c(.5, .75, .9), function(tau) {
        z <- a[a$hour_local == h & a$quantile == tau, ]
        sprintf("%.3f [%.3f, %.3f], p=%.3f%s", z$slope_mph_decade, z$ci_low,
                z$ci_high, z$p_bh_all, if (isTRUE(z$significant_05)) " **significant**" else "")
      }, character(1))
      lines <- c(lines, paste0("| ", h, ":00 | ", paste(entries, collapse = " | "), " |"))
    }
  }
  lines <- c(lines, "", "## Method and limitations", "",
    "Quantile regression is fitted directly to individual hourly sustained-wind observations.",
    "The all-date model uses a linear year term plus annual and semiannual sine/cosine terms, separately by hour and quantile.",
    "The date-window model uses a linear year term without seasonal terms.",
    sprintf("Cluster wild-gradient bootstrap standard errors use %d replicates and calendar-year clusters, with t(%d) tests and pointwise intervals.",
            results$bootstrap_replicates[1], results$year_clusters[1] - 1L),
    "The significance criterion adjusts p-values across all 42 tests; the confidence intervals are not adjusted.",
    "A nonsignificant slope means the data do not establish a linear trend; it does not prove no change occurred.",
    "Year clusters allow within-year dependence but assume independence across years. Rounded wind reports and calm winds produce ties, making inference approximate.",
    "A station trend cannot by itself distinguish climate change from changes in instruments, observing practice, or station exposure.", "",
    "Method reference: [quantreg bootstrap documentation](https://search.r-project.org/CRAN/refmans/quantreg/html/boot.rq.html).", "",
    "Full numerical results: [quantile_trends.csv](quantile_trends.csv).")
  writeLines(lines, file)
}

plot_tlh_wind_trends <- function(results, file) {
  grDevices::png(file, width = 1800, height = 850, res = 160)
  on.exit(grDevices::dev.off())
  old <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1), oma = c(1, 0, 4, 0))
  on.exit(par(old), add = TRUE, after = FALSE)
  limits <- range(c(results$ci_low, results$ci_high, 0), finite = TRUE)
  for (tau in c(.5, .75, .9)) {
    a <- results[results$quantile == tau, ]
    color <- ifelse(a$significant_05, "#177C67", "#555E6B")
    plot(a$hour_local, a$slope_mph_decade, ylim = limits, xaxt = "n", pch = 19,
         col = color, xlab = "Local hour (Eastern)", ylab = "Trend (mph per decade)",
         main = paste0(100 * tau, "th percentile"))
    axis(1, at = 9:15, labels = c("9a", "10a", "11a", "12p", "1p", "2p", "3p"))
    abline(h = 0, lty = 2, col = "#AAAAAA")
    arrows(a$hour_local, a$ci_low, a$hour_local, a$ci_high,
           angle = 90, code = 3, length = .05, col = color)
  }
  label <- if (results$scope[1] == "all_dates") "All dates, adjusted for seasonality" else
    paste0(results$month_day[1], " ±", results$window_days[1], " days")
  mtext(paste0("Tallahassee quantile regression | ", label, " | ",
               results$start_year[1], "–", results$end_year[1]), outer = TRUE, line = 2, cex = 1.2)
  mtext("Bars: pointwise 95% intervals with year-cluster bootstrap SE; green: BH-adjusted p < 0.05 (42 tests)",
        outer = TRUE, line = .5, cex = .8)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  day <- if (length(args)) args[1] else "10-24"
  B <- if (length(args) >= 2L) as.integer(args[2]) else 999L
  result <- tlh_wind_trends(day = day, B = B)
  print(result[c("scope", "hour_local", "quantile", "slope_mph_decade", "ci_low",
                 "ci_high", "p_value", "p_bh_all", "significant_05")], row.names = FALSE)
}
