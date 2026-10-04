# Due to lack of monthly telco data availability clusters 1-4 from Birbilas (2025) 
# were expanded using simulations defined bellow. Custers 5-11 were simulated using 
# cluster 2 data. Details: 
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  | Cluster │ Model   │                Design                │
#  │         │  from   │                                      │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl1     │ cl1     │ no change                            │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl2     │ cl2     │ sqrt(t) jump +25%, months 37–46      │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl3     │ cl3     │ sqrt(t) jump −25%, months 37–49,     │
#  │         │         │ plus outlier months and days         │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl4     │ cl4     │ sqrt(t) jump from month 37,          │
#  │         │         │ deepening 2.5% per month             │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl5     │ cl2     │ rising trend ×1.00 → ×1.30           │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl6     │ cl2     │ falling trend ×1.00 → ×0.70          │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl7     │ cl2     │ higher level ×1.50                   │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl8     │ cl2     │ lower level ×0.60                    │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl9     │ cl2     │ 5% outlier days at 2–4× plus 4       │
#  │         │         │ outlier months at ×1.5               │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl10    │ cl2     │ high and falling ×1.50 → ×1.05       │
#  ├─────────┼─────────┼──────────────────────────────────────┤
#  │ cl11    │ cl2     │ low and rising ×0.60 → ×0.84         │
#  ├─────────┼─────────┼──────────────────────────────────────┤

# Build data/telco_far1.rda: simulated telco daily voice-usage curves that
# follow a functional autoregressive FAR(1) process, 11 clusters x 50 months
# x 28 days. Documented in R/data.R (?telco_far1).
#
# The model parameters in FITS (daily profile mu, month-level sd s_level,
# day-level sd s_day) were fitted on cluster-level averages of a private
# telco bundle-consumption sample (4 clusters, 31 bundle periods x 28 days
# of average daily voice minutes). Only these aggregates are needed, so the
# dataset is fully reproducible from this script.
#
# Model (per cluster, month j = 1..50, day d = 1..28, t_d = (d - 1) / 27):
#   Y[j, d] = f[j] * (mu[d] + Z[j, d])
#             + jump[j] * mean(mu) * sqrt(t_d) / mean(sqrt(t)),   clipped at 0
#   Z[j, d] = beta(t_d) * Z[j-1, d] + sqrt(1 - beta(t_d)^2) * (a[j] + e[j, d])
#   beta(t) = 0.7 * sqrt(t)  (far1_betas$beta1),
#   a[j] ~ N(0, s_level^2), e[j, d] ~ N(0, s_day^2), 100 burn-in months.
#
# Usage: Rscript data-raw/telco_far1.R   (from the package root)

# ---- Design ----------------------------------------------------------------
SEED       <- 2026L
N_MONTHS   <- 50L
FAR_BETA   <- function(t) 0.7 * sqrt(t)
BURN_IN    <- 100L
CL2_WINDOW <- 37:46;  CL2_JUMP <- 0.25     # sqrt(t) jump, month-average size
CL3_WINDOW <- 37:49;  CL3_JUMP <- -0.25
CL3_N_OUTLIERS <- 3L; CL3_OUTLIER_X <- 3   # outlier days: 3 x normal usage
CL3_N_OUT_MONTHS <- 3L; CL3_OUT_MONTH_FACTOR <- 1.5
CL4_START  <- 37L;    CL4_SLOPE <- 0.025   # jump deepens 2.5% per month
NEW_BASE   <- "cl2"                        # model for cl5-cl11
TREND_UP   <- c(1.00, 1.30)                # cl5
TREND_DOWN <- c(1.00, 0.70)                # cl6
LEVEL_HIGH <- 1.50                         # cl7, cl10 start
LEVEL_LOW  <- 0.60                         # cl8, cl11 start
CL9_DAY_SHARE <- 0.05; CL9_DAY_X <- c(2, 4)
CL9_N_OUT_MONTHS <- 4L; CL9_OUT_MONTH_FACTOR <- 1.5
CL10_TREND <- c(LEVEL_HIGH, LEVEL_HIGH * TREND_DOWN[2])
CL11_TREND <- c(LEVEL_LOW,  LEVEL_LOW  * 1.40)

# ---- Fitted parameters -----------------------------------------------------
FITS <- list(cl1 = list(mu = c(8.8865384085383887, 8.6656685540823055, 
8.7559885991082478, 8.7563726120019645, 8.9450219533721587, 8.8335826682951826, 
8.7929218897150747, 9.0824943197904826, 8.8614148686037204, 8.8143335731041681, 
8.7687948629578045, 8.8398342687764568, 8.8732443645114643, 8.6940220850286991, 
8.8097501857346288, 8.6505761735906503, 8.7217885235144266, 8.7289671123725725, 
8.681397256280512, 8.6703180682253773, 8.7633891356335809, 8.7151790788389931, 
8.8702022284583961, 8.6545281494233208, 8.9598701984673212, 8.9280326664466259, 
8.8135081033049474, 8.689082039219393), s_level = 0.12319048852578349, 
    s_day = 0.33093149100900693), cl2 = list(mu = c(9.8807838363751106, 
10.296973586140815, 10.411896246745147, 10.594033404458177, 10.518158358908632, 
10.516315645794954, 10.538491250984682, 10.517065411048764, 10.506889430772707, 
10.583191017655382, 10.505785297299141, 10.524001582557457, 10.572264275273849, 
10.55693534899885, 10.495532279041054, 10.501831955282517, 10.517286953249153, 
10.496516786576773, 10.525874788335601, 10.566663590065781, 10.529175726166866, 
10.468920339221226, 10.427721544927861, 10.464320273920823, 10.39901945340617, 
10.451451638819524, 10.43760966641762, 9.9471244282271023), s_level = 0.95062128681486568, 
    s_day = 0.39553430392436978), cl3 = list(mu = c(8.1090702269059136, 
8.4077718747136387, 8.4671836085641505, 8.5837153401649715, 8.6769464393670948, 
8.7737447793951784, 8.7562779105058191, 8.7803056159668884, 8.7751514479414503, 
8.7293829610926288, 8.7427541398893105, 8.7394000020313882, 8.7456052047675801, 
8.7218278616140257, 8.7624409207120895, 8.7061976963999932, 8.7392460275956605, 
8.735767086416887, 8.6749977022614821, 8.6974029110927873, 8.6588923667258744, 
8.7149499324092652, 8.7599726269277642, 8.6263300951326691, 8.7463161178708866, 
8.7004600785796899, 8.7061112625906709, 8.094381014821499), s_level = 0.59921903322215841, 
    s_day = 0.40216799851258217), cl4 = list(mu = c(9.3521375292788136, 
9.499975368429304, 9.7982407325489298, 9.9785183726624087, 9.9977601404147407, 
10.006761027470358, 9.8361694689047976, 9.7304636099841844, 9.9136021228104241, 
9.983787031516222, 9.9556175872179793, 9.9647255725767909, 10.032041409407624, 
9.961197285842367, 9.9950969508905771, 9.9417908997728262, 9.8576952651202845, 
9.8540390489433349, 10.035185269393999, 9.9732990387899125, 9.9324448864952544, 
10.03799293788548, 9.9935984204438881, 9.9582901950640412, 9.9894420439013043, 
9.8067400911210854, 9.9939091458384457, 9.541227429291542), s_level = 0.76655552712983022, 
    s_day = 0.41073961984049301))

# ---- Model -----------------------------------------------------------------
simulate_cluster <- function(fit, f, jump = rep(0, length(f))) {
  n <- length(f); D <- length(fit$mu)
  t <- seq(0, 1, length.out = D)
  b <- FAR_BETA(t)
  s <- sqrt(1 - b^2)
  Z <- matrix(0, n, D)
  z <- numeric(D)
  for (j in seq_len(BURN_IN + n)) {
    eps <- s * (rnorm(1, 0, fit$s_level) + rnorm(D, 0, fit$s_day))
    z <- b * z + eps
    if (j > BURN_IN) Z[j - BURN_IN, ] <- z
  }
  Y <- f * sweep(Z, 2, fit$mu, "+")
  shape <- sqrt(t) / mean(sqrt(t))
  Y <- Y + outer(jump * mean(fit$mu), shape)
  Y <- pmax(Y, 0)
  dimnames(Y) <- list(as.character(seq_len(n)), as.character(seq_len(D)))
  Y
}

no_outliers <- function() {
  matrix(integer(0), 0, 2, dimnames = list(NULL, c("month", "day")))
}
pick_days <- function(months, n_days, n) {
  cells <- sample(length(months) * n_days, n)
  cbind(month = months[(cells - 1L) %% length(months) + 1L],
        day   = (cells - 1L) %/% length(months) + 1L)
}
finish <- function(Y, change_months, f, jump, out, out_months) {
  attr(Y, "change_months")  <- change_months
  attr(Y, "multiplier")     <- f
  attr(Y, "jump")           <- jump
  attr(Y, "outliers")       <- out
  attr(Y, "outlier_months") <- out_months
  Y
}

set.seed(SEED)
telco_far1 <- list()

# ---- cl1-cl4: each on its own model, sqrt(t) jumps in cl2-cl4 --------------
zero <- rep(0, N_MONTHS)
jump <- list(
  cl1 = zero,
  cl2 = replace(zero, CL2_WINDOW, CL2_JUMP),
  cl3 = replace(zero, CL3_WINDOW, CL3_JUMP),
  cl4 = c(rep(0, CL4_START - 1L),
          -CL4_SLOPE * (seq_len(N_MONTHS - CL4_START + 1L) - 1L)))
mult <- lapply(jump, function(x) rep(1, N_MONTHS))
window <- list(cl1 = integer(0), cl2 = CL2_WINDOW, cl3 = CL3_WINDOW,
               cl4 = CL4_START:N_MONTHS)
out_months <- list(cl1 = integer(0), cl2 = integer(0), cl4 = integer(0),
                   cl3 = sort(sample(CL3_WINDOW, CL3_N_OUT_MONTHS)))
mult$cl3[out_months$cl3] <- CL3_OUT_MONTH_FACTOR   # outlier month, no jump
jump$cl3[out_months$cl3] <- 0

for (nm in names(mult)) {
  fit <- FITS[[nm]]
  Y <- simulate_cluster(fit, mult[[nm]], jump[[nm]])
  out <- no_outliers()
  if (nm == "cl3") {
    out <- pick_days(setdiff(CL3_WINDOW, out_months$cl3), ncol(Y),
                     CL3_N_OUTLIERS)
    Y[out] <- CL3_OUTLIER_X * fit$mu[out[, "day"]]
  }
  telco_far1[[nm]] <- finish(Y, window[[nm]], mult[[nm]], jump[[nm]], out,
                             out_months[[nm]])
}

# ---- cl5-cl11: no jump, cl2's model, level multiplier only -----------------
trend <- function(from_to) seq(from_to[1], from_to[2], length.out = N_MONTHS)
mult_new <- list(cl5  = trend(TREND_UP),
                 cl6  = trend(TREND_DOWN),
                 cl7  = rep(LEVEL_HIGH, N_MONTHS),
                 cl8  = rep(LEVEL_LOW, N_MONTHS),
                 cl9  = rep(1, N_MONTHS),
                 cl10 = trend(CL10_TREND),
                 cl11 = trend(CL11_TREND))
base_fit <- FITS[[NEW_BASE]]
for (nm in names(mult_new)) {
  f <- mult_new[[nm]]
  om <- integer(0)
  if (nm == "cl9") {
    repeat {                       # isolated outlier months, never adjacent
      om <- sort(sample(N_MONTHS, CL9_N_OUT_MONTHS))
      if (all(diff(om) > 1L)) break
    }
    f[om] <- CL9_OUT_MONTH_FACTOR
  }
  Y <- simulate_cluster(base_fit, f)
  out <- no_outliers()
  if (nm == "cl9") {
    day_months <- setdiff(seq_len(N_MONTHS), om)
    n_out <- round(CL9_DAY_SHARE * length(day_months) * ncol(Y))
    out <- pick_days(day_months, ncol(Y), n_out)
    Y[out] <- runif(n_out, CL9_DAY_X[1], CL9_DAY_X[2]) *
      base_fit$mu[out[, "day"]]
  }
  telco_far1[[nm]] <- finish(Y, integer(0), f, rep(0, N_MONTHS), out, om)
}

dir.create("data", showWarnings = FALSE)
save(telco_far1, file = file.path("data", "telco_far1.rda"), compress = "xz",
     version = 2)
