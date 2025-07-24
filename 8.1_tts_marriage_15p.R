#---- Sensitivity analysis - cohabitations from age 15 ----
#---- Transitions from cohabitation to marriage ----

#---- Load the necessary libraries and functions ----

library(data.table)
library(TwoTimeScales)

# ---- Load the data ----
load("one_cohabitation15.Rda")

## Separate by gender and east/west
WomWest <- one_cohab_15p[sex == 2 & east == 0]
WomEast <- one_cohab_15p[sex == 2 & east == 1]
MenWest <- one_cohab_15p[sex == 1 & east == 0]
MenEast <- one_cohab_15p[sex == 1 & east == 1]


# ---- Two time scales analysis without covariates ----

#* Men east ----
du <- ds <- .5 # half-year
data2d <- prepare_data(
  u = MenEast$agecoh_y,
  s_in = MenEast$dur_entry,
  s_out = MenEast$durcoh_y,
  events = as.numeric(MenEast$event_type == "marriage"),
  ds = ds,
  du = du
)
range(MenEast$agecoh_y)
range(MenEast$durcoh_y)

optimal_model <- fit2ts(data2d,
                        Bbases_spec = list(
                          nseg_u = 14,
                          nseg_s = 12,
                          min_u = 15, max_u = 47,
                          min_s = 0, max_s = 30
                        )
)

# let's add the intervals, the B-splines and their specification
# to the saved results
r_mar_nocov_15p_me <- list(
  mod_no_cov = optimal_model,
  data2d = data2d
)
save(r_mar_nocov_15p_me, file = "R_marriage_nocov_15plus_me.Rda")

#* Men West ----
du <- ds <- .5 # half-year
data2d <- prepare_data(
  u = MenWest$agecoh_y,
  s_in = MenWest$dur_entry,
  s_out = MenWest$durcoh_y,
  events = as.numeric(MenWest$event_type == "marriage"),
  ds = ds,
  du = du
)
range(MenWest$agecoh_y)
range(MenWest$durcoh_y)
optimal_model <- fit2ts(data2d,
                        Bbases_spec = list(
                          nseg_u = 14,
                          nseg_s = 9,
                          min_u = 15, max_u = 47,
                          min_s = 0, max_s = 25
                        )
)

# let's add the intervals, the B-splines and their specification
# to the saved results
r_mar_nocov_15p_mw <- list(
  mod_no_cov = optimal_model,
  data2d = data2d
)
save(r_mar_nocov_15p_mw, file = "R_marriage_nocov_15plus_mw.Rda")

#* Women East ----
du <- ds <- .5 # half-year
data2d <- prepare_data(
  u = WomEast$agecoh_y,
  s_in = WomEast$dur_entry,
  s_out = WomEast$durcoh_y,
  events = as.numeric(WomEast$event_type == "marriage"),
  ds = ds,
  du = du
)
range(WomEast$agecoh_y)
range(WomEast$durcoh_y)
optimal_model <- fit2ts(data2d,
                        Bbases_spec = list(
                          nseg_u = 14,
                          nseg_s = 9,
                          min_u = 15, max_u = 47,
                          min_s = 0, max_s = 25
                        )
)

# let's add the intervals, the B-splines and their specification
# to the saved results
r_mar_nocov_15p_we <- list(
  mod_no_cov = optimal_model,
  data2d = data2d
)
save(r_mar_nocov_15p_we, file = "R_marriage_nocov_15plus_we.Rda")

#* Women West ----
du <- ds <- .5 # half-year
data2d <- prepare_data(
  u = WomWest$agecoh_y,
  s_in = WomWest$dur_entry,
  s_out = WomWest$durcoh_y,
  events = as.numeric(WomWest$event_type == "marriage"),
  ds = ds,
  du = du
)
range(WomWest$agecoh_y)
range(WomWest$durcoh_y)
optimal_model <- fit2ts(data2d,
                        Bbases_spec = list(
                          nseg_u = 14,
                          nseg_s = 10,
                          min_u = 15, max_u = 47,
                          min_s = 0, max_s = 27
                        )
)

# let's add the intervals, the B-splines and their specification
# to the saved results
r_mar_nocov_15p_ww <- list(
  mod_no_cov = optimal_model,
  data2d = data2d
)
save(r_mar_nocov_15p_ww, file = "R_marriage_nocov_15plus_ww.Rda")
