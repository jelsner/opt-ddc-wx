# Run from the project root: Rscript scripts/test_tallahassee_wind_exceedance.R
source("scripts/tallahassee_wind_exceedance.R")

stopifnot(identical(tlh_years("2026-10-24"), 1996:2025),
          tlh_day(60) == "03-01", tlh_day(60, 2000) == "02-29",
          tlh_day("2026-07-04") == "07-04", tlh_day("02-29") == "02-29")
stopifnot(inherits(try(tlh_day("02-30"), silent = TRUE), "try-error"))

# Same local 09:00 in winter and summer, duplicate valid/missing reports,
# invalid speeds, and reports outside the requested hours.
raw <- data.frame(station = "TLH", valid = c(
  "2025-01-01 13:53", "2025-07-01 12:53", "2025-07-01 12:54",
  "2025-07-01 13:53", "2025-07-01 14:53", "2025-07-01 11:53",
  "2025-07-01 19:53"), sknt = c(0, 10, NA, 20, -1, 90, 90))
p <- tlh_prepare(raw, 1996:2025)
stopifnot(nrow(p) == 4L, all(p$hour == c(9, 9, 10, 11)),
          p$wind_knots[2] == 10, is.na(p$wind_knots[4]))
r <- suppressWarnings(tlh_wind_exceedance("07-01", c(0, 10, 20), "knots",
  window_days = 0, end_year = 2025, observations = raw))
stopifnot(identical(r$probability[r$hour_local == 9], c(0, 1, 1)),
          all(is.na(r$probability[r$hour_local == 11])),
          all(r$expected_n == 30L), all(r$n[r$hour_local == 9] == 1L))

# Five dates around Jan 1 wrap around New Year, with one vote per year/date/hour.
dates <- c("2024-12-30", "2024-12-31", "2025-01-01", "2025-01-02", "2025-01-03")
raw2 <- data.frame(station = "TLH", valid = paste(dates, "13:53"),
                   sknt = c(0, 5, 10, 15, 20))
r <- suppressWarnings(tlh_wind_exceedance("01-01", c(0, 10, 20), "knots",
  end_year = 2025, observations = raw2))
stopifnot(all(r$expected_n == 150L),
          isTRUE(all.equal(r$probability[r$hour_local == 9], c(0.2, 0.6, 1))),
          all(r$probability + r$exceedance_probability == 1, na.rm = TRUE),
          all(r$n[r$hour_local == 9] == 5L))

# Leap-day handling: exact Feb 29 only has eight possible years in 1996–2025.
# Feb 28 ±2 days has five actual dates in each year, including leap years.
r <- suppressWarnings(tlh_wind_exceedance("02-29", 5, window_days = 0,
  end_year = 2025, observations = raw2))
stopifnot(all(r$expected_n == 8L))
r <- suppressWarnings(tlh_wind_exceedance("02-28", 5,
  end_year = 2025, observations = raw2))
stopifnot(all(r$expected_n == 150L))
r <- suppressWarnings(tlh_wind_exceedance("02-29", 5,
  end_year = 2025, observations = raw2))
stopifnot(all(r$expected_n == 40L))

# The mph table includes a threshold above the largest selected wind.
r <- suppressWarnings(tlh_wind_exceedance("01-01", end_year = 2025,
  observations = transform(raw2, sknt = 40)))
stopifnot(max(r$threshold) == ceiling(40 * 1852 / 1609.344),
          all(diff(r$probability[r$hour_local == 9]) >= 0),
          tail(r$probability[r$hour_local == 9], 1) == 1)
cat("All Tallahassee wind tests passed.\n")

source("scripts/tallahassee_wind_quantiles.R")
q <- suppressWarnings(tlh_wind_quantiles("01-01", units = "knots", end_year = 2025,
                                        observations = raw2, plot = FALSE))
stopifnot(identical(q$wind_speed[q$hour_local == 9], c(20, 15, 10, 5, 0)),
          all(q$n[q$hour_local == 9] == 5L), all(q$expected_n == 150L),
          all(is.na(q$wind_speed[q$hour_local == 10])))
mph <- suppressWarnings(tlh_wind_quantiles("01-01", probs = .5, end_year = 2025,
                                          observations = raw2, plot = FALSE))
stopifnot(isTRUE(all.equal(mph$wind_speed[1], 10 * 1852 / 1609.344)))
cat("All Tallahassee quantile tests passed.\n")
