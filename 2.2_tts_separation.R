#---- Transitions from cohabitation to separation ----

#---- Load the necessary libraries and functions ----

library(data.table)
library(TwoTimeScales)

# ---- Load the data ----
load("one_cohabitation18.Rda")

# ---- Separate by gender and east/west
WomWest <- one_cohab[sex == 2 & east == 0]
WomEast <- one_cohab[sex == 2 & east == 1]
MenWest <- one_cohab[sex == 1 & east == 0]
MenEast <- one_cohab[sex == 1 & east == 1]

# ---- Two time scales analysis without covariates ----

#* Men east ----
du <- ds <- .5 # half-year

range(MenEast$agecoh_y)
range(MenEast$durcoh_y)
data2d <- prepare_data(
  u = MenEast$agecoh_y,
  s_in = MenEast$dur_entry,
  s_out = MenEast$durcoh_y,
  events = as.numeric(MenEast$event_type == "separation"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47
)
optimal_model <- fit2ts(data2d,
  Bbases_spec = list(
    nseg_u = 11,
    nseg_s = 10,
    min_u = 18, max_u = 47,
    min_s = 0, max_s = 26
  )
)

# let's add the intervals, the B-splines and their specification
# to the saved results
r_sep_nocov_18p_me <- list(
  mod_no_cov = optimal_model,
  data2d = data2d
)
save(r_sep_nocov_18p_me, file = "R_int2d_separation_me.Rda")

#* Men West ----
du <- ds <- .5 # half-year
range(MenWest$agecoh_y)
range(MenWest$durcoh_y)
data2d <- prepare_data(
  u = MenWest$agecoh_y,
  s_in = MenWest$dur_entry,
  s_out = MenWest$durcoh_y,
  events = as.numeric(MenWest$event_type == "separation"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47
)
optimal_model <- fit2ts(data2d,
  Bbases_spec = list(
    nseg_u = 11,
    nseg_s = 9,
    min_u = 18, max_u = 47,
    min_s = 0, max_s = 25
  )
)

# let's add the intervals, the B-splines and their specification
# to the saved results
r_sep_nocov_18p_mw <- list(
  mod_no_cov = optimal_model,
  data2d = data2d
)
save(r_sep_nocov_18p_mw, file = "R_int2d_separation_mw.Rda")

#* Women East ----
du <- ds <- .5 # half-year
range(WomEast$agecoh_y)
range(WomEast$durcoh_y)
data2d <- prepare_data(
  u = WomEast$agecoh_y,
  s_in = WomEast$dur_entry,
  s_out = WomEast$durcoh_y,
  events = as.numeric(WomEast$event_type == "separation"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47
)
optimal_model <- fit2ts(data2d,
  Bbases_spec = list(
    nseg_u = 11,
    nseg_s = 9,
    min_u = 18, max_u = 47,
    min_s = 0, max_s = 25
  )
)

# let's add the intervals, the B-splines and their specification
# to the saved results
r_sep_nocov_18p_we <- list(
  mod_no_cov = optimal_model,
  data2d = data2d
)
save(r_sep_nocov_18p_we, file = "R_int2d_separation_we.Rda")

#* Women West ----
du <- ds <- .5 # half-year
range(WomWest$agecoh_y)
range(WomWest$durcoh_y)
data2d <- prepare_data(
  u = WomWest$agecoh_y,
  s_in = WomWest$dur_entry,
  s_out = WomWest$durcoh_y,
  events = as.numeric(WomWest$event_type == "separation"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47
)

optimal_model <- fit2ts(data2d,
  Bbases_spec = list(
    nseg_u = 11,
    nseg_s = 10,
    min_u = 18, max_u = 47,
    min_s = 0, max_s = 27
  )
)

# let's add the intervals, the B-splines and their specification
# to the saved results
r_sep_nocov_18p_ww <- list(
  mod_no_cov = optimal_model,
  data2d = data2d
)
save(r_sep_nocov_18p_ww, file = "R_int2d_separation_ww.Rda")
