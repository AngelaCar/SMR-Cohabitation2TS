# ---- Transitions from cohabitation to separation ----
# Varying-coefficient models

#---- Load the necessary libraries and functions ----

library(data.table)
library(TwoTimeScales)
library(JOPS)
source("VCM_2TS_fixed.R")

# ---- Load the data ----
load("one_cohabitation18.Rda")

## Separate by gender and east/west
WomWest <- one_cohab[sex == 2 & east == 0]
WomEast <- one_cohab[sex == 2 & east == 1]
MenWest <- one_cohab[sex == 1 & east == 0]
MenEast <- one_cohab[sex == 1 & east == 1]

#* Men east ----
du <- ds <- .5 # half-year
data2d_me <- prepare_data(
  u = MenEast$agecoh_y,
  s_in = MenEast$dur_entry,
  s_out = MenEast$durcoh_y,
  events = as.numeric(MenEast$event_type == "separation"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47
)

vcm_sep_me <- VCM_2TS_fixed(
  ExEv = data2d_me$bindata,
  binlist = data2d_me$bins,
  nseg = 10,
  bdeg = 3, pord = 2, kappa = 1e-10,
  lambdas = c(10^3.8, 10^5.5)
)

r_sep_vcm_me <- list(
  vcm_sep_me = vcm_sep_me,
  data2d = data2d_me
)
save(r_sep_vcm_me, file = "R_separation_vcm_me.Rda")

#* Men west ----
du <- ds <- .5 # half-year
data2d_mw <- prepare_data(
  u = MenWest$agecoh_y,
  s_in = MenWest$dur_entry,
  s_out = MenWest$durcoh_y,
  events = as.numeric(MenWest$event_type == "separation"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47,
  min_s = 0, max_s = 25
)

vcm_sep_mw <- VCM_2TS_fixed(
  ExEv = data2d_mw$bindata,
  binlist = data2d_mw$bins,
  nseg = 9, bdeg = 3, pord = 2, kappa = 1e-10,
  lambdas = c(10^4, 10^6.5)
)
r_sep_vcm_mw <- list(
  vcm_sep_mw = vcm_sep_mw,
  data2d = data2d_mw
)
save(r_sep_vcm_mw, file = "R_separation_vcm_mw.Rda")

#* Women east ----
du <- ds <- .5 # half-year
data2d_we <- prepare_data(
  u = WomEast$agecoh_y,
  s_in = WomEast$dur_entry,
  s_out = WomEast$durcoh_y,
  events = as.numeric(WomEast$event_type == "separation"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47
)

vcm_sep_we <- VCM_2TS_fixed(
  ExEv = data2d_we$bindata,
  binlist = data2d_we$bins,
  nseg = 9,
  bdeg = 3, pord = 2, kappa = 1e-10,
  lambdas = c(10^3.3, 10^5.3)
)

r_sep_vcm_we <- list(
  vcm_sep_we = vcm_sep_we,
  data2d = data2d_we
)
save(r_sep_vcm_we, file = "R_separation_vcm_we.Rda")

#* Women West ----
du <- ds <- .5 # half-year
data2d_ww <- prepare_data(
  u = WomWest$agecoh_y,
  s_in = WomWest$dur_entry,
  s_out = WomWest$durcoh_y,
  events = as.numeric(WomWest$event_type == "separation"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47,
  min_s = 0, max_s = 27
)

vcm_sep_ww <- VCM_2TS_fixed(
  ExEv = data2d_ww$bindata,
  binlist = data2d_ww$bins,
  nseg = 10, bdeg = 3, pord = 2, kappa = 1e-10,
  lambdas = c(10^5.5, 10^5)
)

r_sep_vcm_ww <- list(
  vcm_sep_ww = vcm_sep_ww,
  data2d = data2d_ww
)
save(r_sep_vcm_ww, file = "R_separation_vcm_ww.Rda")
