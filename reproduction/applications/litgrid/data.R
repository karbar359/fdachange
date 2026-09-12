#' Litgrid (Lithuanian electricity-balancing market) data loader.
#'
#' Reads raw CSV exports from `<raw_dir>/litgrid/` (downloaded from
#' https://baltic.transparency-dashboard.eu/) and reshapes them into a
#' functional-data matrix whose columns are week-long curves sampled at
#' hour-of-week resolution (168 hours/week).
#'
#' Pipeline (matching `code/practical_applications/litgrid/litgrid_hourly.R`):
#'   1. Concatenate all `balancing-energy-prices.*.csv` files.
#'   2. Filter rows with area == "Lithuania".
#'   3. Derive value = Upward - Downward (decimal comma -> numeric).
#'   4. Aggregate to hourly by summing same-hour rows.
#'   5. Interpolate missing hours, complete weeks, merge year-boundary weeks.
#'   6. Pivot to wide: rows = year_week, cols = 1..168 (hour_of_week).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(zoo)
})

#' Load and reshape paper 1's primary Litgrid dataset (balancing prices).
#'
#' @param raw_dir path to the `<repo>/.raw/litgrid/` directory.
#' @param dedup_exports if TRUE (default) take one value per hour when the raw
#'   exports overlap. The 2018 export runs to 2019-10-01 and the 2019 export
#'   starts at 2019-01-01, so ~6,551 hours in Jan-Sep 2019 are present twice;
#'   the original archive script summed them, which silently doubled those
#'   hours. Pass FALSE to reproduce the printed dissertation tables, which
#'   inherit that behaviour (see `reproduce_tables.R`).
#' @param fix_week_labels if TRUE (default) key each Monday-start week by its
#'   ISO year-week (via the week's Thursday) and order rows chronologically.
#'   The legacy labelling (`year = min(year), week = min(isoweek)` after
#'   grouping) mislabels boundary-spanning W01 weeks: the 2019-12-30 week is
#'   labelled 201901 and string-sorts to position 52, a genuine break in the
#'   time order the scan statistics assume, and the true 2018-12-31 week
#'   collides with the truncated 2018-01-01 stub into a 335-hour key that the
#'   168-hour filter then silently drops. FALSE reproduces the print.
#' @param stamp "from" (default) keys each settlement hour by its interval
#'   START, so column 1 of a weekly curve is Monday 00:00-01:00. The legacy
#'   convention ("to") keys by the interval END, which rotates every weekly
#'   curve one hour: column 1 becomes the LAST hour of the preceding week.
#' @return numeric matrix with rows = year_week label, cols = hour_of_week (168).
load_litgrid_balancing_prices <- function(raw_dir, dedup_exports = TRUE,
                                          fix_week_labels = TRUE,
                                          stamp = c("from", "to")) {
  stamp <- match.arg(stamp)
  files <- list.files(raw_dir,
                      pattern = "^balancing-energy-prices.*\\.csv$",
                      full.names = TRUE)
  if (!length(files)) {
    stop("No balancing-energy-prices CSVs under ", raw_dir)
  }
  raw <- read_delim(files, delim = ";", escape_double = FALSE,
                    trim_ws = TRUE, id = "id", show_col_types = FALSE)
  d <- subset(raw, area == "Lithuania")
  d$Upward.n   <- as.numeric(gsub(",", ".", d$Upward))
  d$Downward.n <- as.numeric(gsub(",", ".", d$Downward))
  d$value      <- d$Upward.n - d$Downward.n
  d$datetime   <- if (stamp == "from") ymd_hm(d$datetime_from) else ymd_hm(d$datetime_to)
  d_hourly <- d %>%
    mutate(datetime_hour = floor_date(datetime, unit = "hour"))
  if (dedup_exports) {
    # Collapse repeats within one export first (the DST fall-back hour appears
    # twice under a local-time stamp), preserving the original behaviour there,
    # then take one value per hour across the overlapping exports instead of
    # adding them together. Only 3 of the 6,557 duplicated hours disagree
    # between exports, so the mean is a rounding-level choice.
    d_hourly <- d_hourly %>%
      group_by(id, datetime_hour) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      group_by(datetime_hour) %>%
      summarise(value = mean(value), .groups = "drop")
  } else {
    d_hourly <- d_hourly %>%
      group_by(datetime_hour) %>%
      summarise(value = sum(value), .groups = "drop")
  }
  reshape_hourly_to_week_matrix(d_hourly, "datetime_hour", "value",
                                fix_week_labels = fix_week_labels)
}

#' Convert a long (datetime, value) data frame into a (year_week x 168) matrix.
#' With `fix_week_labels = FALSE` this ports
#' `helper_functions.R::prep.data.into.matrix` verbatim, defects included.
reshape_hourly_to_week_matrix <- function(data, datetime_column, value_column,
                                          fix_week_labels = FALSE) {
  datetime_col <- sym(datetime_column)
  value_col    <- sym(value_column)

  df_complete <- data %>%
    tidyr::complete(
      !!datetime_col := seq(min(!!datetime_col), max(!!datetime_col),
                            by = "hour")
    ) %>%
    mutate(missing_flag = ifelse(is.na(!!value_col), 1, 0)) %>%
    mutate(!!value_col := zoo::na.approx(!!value_col, na.rm = FALSE)) %>%
    mutate(value = !!value_col)

  df_dt <- df_complete %>%
    mutate(
      year         = year(!!datetime_col),
      week_of_year = sprintf("%02d", isoweek(!!datetime_col)),
      day_of_week  = wday(!!datetime_col, week_start = 1),
      hour_of_day  = hour(!!datetime_col),
      week_start   = floor_date(!!datetime_col, unit = "week", week_start = 1),
      week_end     = ceiling_date(!!datetime_col, unit = "week",
                                  week_start = 1) - days(1)
    )

  if (fix_week_labels) {
    # Key each Monday-start week by its ISO year-week, taken from the week's
    # Thursday (the ISO-8601 rule), and order rows by week_start. Labels are
    # then bijective with calendar weeks: no boundary-week teleport, no
    # 201801-style collision.
    df_merged_weeks <- df_dt %>%
      mutate(hour_of_week = (day_of_week - 1) * 24 + hour_of_day + 1) %>%
      mutate(year_week = sprintf("%d%02d",
                                 isoyear(week_start + days(3)),
                                 isoweek(week_start + days(3))))

    week_counts <- df_merged_weeks %>%
      count(year_week, name = "n_hours") %>%
      filter(n_hours == 168)

    df_full_weeks <- df_merged_weeks %>%
      semi_join(week_counts, by = "year_week")

    ws_by_label <- df_full_weeks %>%
      distinct(year_week, week_start) %>%
      arrange(week_start)
    if (anyDuplicated(ws_by_label$year_week)) {
      stop("ISO week labels are not unique; labelling invariant broken.")
    }

    df_wide <- df_full_weeks[, c("year_week", "hour_of_week", "value")] %>%
      pivot_wider(names_from = hour_of_week, values_from = value)

    m <- as.data.frame(df_wide)
    rownames(m) <- m$year_week
    m <- as.matrix(m[, -1L])
    return(m[ws_by_label$year_week, , drop = FALSE])
  }

  df_merged_weeks <- df_dt %>%
    group_by(week_start, week_end) %>%
    mutate(week_of_year = min(week_of_year), year = min(year)) %>%
    ungroup() %>%
    mutate(hour_of_week = (day_of_week - 1) * 24 + hour_of_day + 1) %>%
    unite("year_week", year, week_of_year, sep = "", remove = FALSE)

  week_counts <- df_merged_weeks %>%
    group_by(year, week_of_year) %>%
    summarise(n_hours = n(), .groups = "drop") %>%
    filter(n_hours == 168)

  df_full_weeks <- df_merged_weeks %>%
    inner_join(week_counts, by = c("year", "week_of_year"))

  df_wide <- df_full_weeks[, c("year_week", "hour_of_week", "value")] %>%
    pivot_wider(names_from = hour_of_week, values_from = value) %>%
    arrange(year_week)

  m <- as.data.frame(df_wide)
  rownames(m) <- m$year_week
  m <- as.matrix(m[, -1L])
  m
}
