# ---- Sensitivity Analysis - comparison with Hiekel et al., 2015 ----
# ---- libraries and stata dataset ----
library(readstata13)
library(Hmisc)
library(data.table)
library(TwoTimeScales)

# ---- 1: Data preparation ----
# ----- data cleaning ----- 
data <- read.dta13("anchor_w1-4.dta", # prepared externally in Stata
                   convert.factors = FALSE,
                   replace.strl = TRUE)

# remove individuals with only one observation
data <- as.data.table(data)
data2 <- data[, .SD[.N > 1], by = id]

# select only cohorts 1971-73 and 1981-83 (cohort 2, 3)
data3 <- data2[cohort %in% 2:3]

# select only heterosexual individuals
data4 <- data3[homosex == 0]

# only individuals in a cohabitation at first interview (no previous marriage)
coh_w1 <- data4[wave == 1 & relstat %in% c(3,8,11)]
data5 <- data4[id %in% coh_w1$id]

describe(data5$id)

data6 <- data5[(relstat %in% c(1,3,4,6,8,9,11))]

# Event definition
table(data6[wave == 1, nmar])

data6[, inmar := ifelse(relstat == 4, 1, 0)]
data6[, separated := ifelse(relstat %in% c(1,6,9), 1, 0)]
#View(data6[, c("id", "wave", "relstat", "reldur", "cohabdur", "inmar", "separated")])

# identify only first event
data6[, cumevent_mar := cumsum(inmar), by = id]
data6[, cumevent_sep := cumsum(separated), by = id]
#View(data6[, c("id", "wave", "relstat", "reldur",  "cohabdur", "mardur", "inmar", "cumevent_mar", "cumevent_sep")])

# keep only until no events before first event
data7 <- data6[cumevent_mar <= 1]
data7 <- data7[cumevent_sep <= 1]
setorder(data7, id, wave)
# View(data7[, c("id", "wave", "relstat", "reldur", "cohabdur", "mardur", "inmar",
#                "cumevent_mar", "cumevent_sep")])

# Define start and stop
data7[, datint_next := Hmisc::Lag(dateint, -1), by = "id"]
# View(data7[, c("id", "wave", "relstat", "reldur", "inmar", "dateint", 
#                "datint_next", "cohabdur", "mardur")])
data7[, start := dateint]
data7[inmar == 1, stop := dateint - mardur]
data7[, stop_next := Hmisc::Lag(stop, -1), by = "id"]
# View(data7[, c("id", "wave", "relstat", "reldur", "inmar", "dateint", 
#                "datint_next", "cohabdur", "mardur", "stop_next")])
data7[!is.na(stop_next), stop := stop_next , by = "id"]
data7[, mar_lag := Hmisc::Lag(inmar, -1), by = "id"]
data7[, sep_lag := Hmisc::Lag(separated, -1), by = "id"]

# View(data7[, c("id", "wave", "dateint", "datint_next", 
#                "cohabdur", "inmar", "separated", "stop_next", "stop",
#                "mar_lag", "sep_lag")])
data7[is.na(stop), stop := datint_next]

data8 <- data7[!is.na(stop)]
nrow(data8[start > stop])
# View(data8[start > stop,
#            c("id", "wave", "dateint", "datint_next", 
#              "cohabdur", "inmar", "separated", "stop_next", "stop",
#              "mar_lag", "sep_lag")])
data9 <- data8[!(start > stop)]
nrow(data9[start == stop])

data9[start == stop, stop := stop + 0.5]

# Now translate start and stop in ages
data9[, dob := ((doby_gen - 1900) * 12 + dobm_gen - 1)]

data9[, age_start := (start - dob)/12]
data9[, age_stop := (stop - dob)/12]
# covariates
#data9[marstat == 2, marbeg := dateint - mardur ] 

#data9[marstat == 2, agemar := marbeg - dob]
#data9[marstat == 2, agemar_y := agemar / 12]

# Education
data9[, educ_cat := ifelse(isced2 == -7, "Unknown", NA)]
data9[, educ_cat := ifelse(isced2 == 1 | isced2 == 2 | isced2 == 3,
                           "Low", educ_cat)]
data9[, educ_cat := ifelse(isced2 == 4 | isced2 == 5 | isced2 == 6,
                           "Medium", educ_cat)]
data9[, educ_cat := ifelse(isced2 == 7 | isced2 == 8 ,
                           "High", educ_cat)]
data9[, educ_cat := ifelse(isced2 == 0, "Enrolled", educ_cat)]

describe(data9$educ_cat)
# remove education unknown
data10 <- data9[!educ_cat == "Unknown"]

describe(data10$east)

# employment
describe(data10$lfs)
data10[, empl := ifelse(lfs == -7, "Unknown", NA)]
data10[, empl := ifelse(lfs == 1, "InEducation", empl)]
data10[, empl := ifelse(lfs %in% 2:8, "Unemployment", empl)]
data10[, empl := ifelse(lfs %in% 9:13, "Employment", empl)]
describe(data10$empl)
# remove employment unknown
data11 <- data10[!empl == "Unknown"]

table(data11$educ_cat, data11$empl)

# religiosity
table(data11$sd30, data11$sd31)
data11[wave == 1, rel := ifelse(sd30 %in% 1:6 & sd31 %in% 1:4, "religious", NA)]
data11[wave == 1, rel := ifelse(sd30 == 7 | (sd30 %in% 1:6 & sd31 %in% 5:6), 
                                "non_religious", rel)]
data11[wave == 1 & is.na(rel), rel := "undetermined"]
describe(data11[wave == 1, rel])
data11[, rel := zoo::na.locf(rel), by = id]
describe(data11[, rel])

# remove undetermined
data12 <- data11[!rel == "undetermined"]

# typology of cohabitation
# marriage plans
describe(data12[, pa11])
data12[, plans := ifelse(pa11 %in% 1:2, "yes", NA)]
data12[, plans := ifelse(pa11 %in% 3:5, "no", plans)]
data12[, plans := ifelse(pa11 %in% c(-4,-3,-2,-1), "undetermined", plans)]
describe(data12[wave == 1, plans])

# values regarding marriage - keep as fixed
data12[wave == 1, value := ifelse(val1i2 %in% 1:2, "disagree", NA)]
data12[wave == 1, value := ifelse(val1i2 %in% c(3,-1), "indifferent", value)]
data12[wave == 1, value := ifelse(val1i2 %in% 4:5, "agree", value)]
# data6[wave == 3, value := ifelse(val1i2 %in% 1:2, "disagree", NA)]
# data6[wave == 3, value := ifelse(val1i2 %in% c(3,-1), "indifferent", value)]
# data6[wave == 3, value := ifelse(val1i2 %in% 4:5, "agree", value)]

data12[(wave == 1 ) & is.na(value), value := "undetermined"]
describe(data12[wave == 1 , value])
# fill in for next waves
data12[, value := zoo::na.locf(value), by = id]

# finally, typology of cohabitation
data13 <- data12[!(wave == 1 & plans == "undetermined")]
data13[wave == 1, coh_type := ifelse(plans == "yes" & value == "agree", "prelude", NA)]
data13[wave == 1, coh_type := ifelse(plans == "yes" & (value == "disagree" | value == "indifferent"), "conformist", coh_type)]
data13[wave == 1, coh_type := ifelse(plans == "no" & value == "agree", "notready", coh_type)]
data13[wave == 1, coh_type := ifelse(plans == "no" & value == "indifferent", "irrelevant", coh_type)]
data13[wave == 1, coh_type := ifelse(plans == "no" & value == "disagree", "reject", coh_type)]
describe(data13[wave == 1, coh_type])
data13[, coh_type := zoo::na.locf(coh_type), by = id]

# kids
describe(data13[, c("nkids", "k1type")])
data13[, kids := ifelse(nkids == 0, "no_child", NA)]
data13[, kids := ifelse(nkids > 0 & 
                          (k1type %in% c(3, 6, 7, 9)), "child_w_partner", kids)]
data13[, kids := ifelse(nkids > 0 & (k1type %in% c(2, 4, 5, 8)),
                        "child_other_partner", kids)]
describe(data13$kids)
table(data13[is.na(kids), nkids])
table(data13[is.na(kids), k1type])

data14 <- data13[!is.na(kids)]

# previous marriage
data15 <- data14[nmar>= 0]
data15[, prevmar := ifelse(nmar > 0, 1, 0)]

# find duration of cohabitation at first interview and at last observation
data15[, observ_ind := 1:.N, by = id]
data15[observ_ind == 1, dur_first_obs := cohabdur]
data15[observ_ind == 1, cohbeg := dateint - cohabdur]
data15[observ_ind == 1, agecoh := cohbeg - dob]
data15[observ_ind == 1, agecoh_y := agecoh/12]
data15[, tot_obs := .N, by = id]
data15[observ_ind == tot_obs, stop_last_obs := stop]

# impute age at last observation for all
data15[, stop_last_obs := zoo::na.locf(stop_last_obs, fromLast = TRUE), by = id]
# if observation span ends with event, this should be recorded for every individual
data15[, marriage := max(mar_lag), by = id]
data15[, separation := max(sep_lag), by = id]
data_baseline <- data15[observ_ind == 1] 
data_baseline[, dur_last_obs := (stop_last_obs - start) + dur_first_obs]

save(data_baseline, file = "data_wide_format-asHiekel2015.Rda")


# ---- 2: Analysis ----
load("data_wide_format-asHiekel2015.Rda")

data_baseline[, censoring_ind := ifelse(marriage == 1, "marriage", NA)]
data_baseline[, censoring_ind := ifelse(separation == 1, "separation", censoring_ind)]
data_baseline[, censoring_ind := ifelse(marriage == 0 & separation == 0, "censored", censoring_ind)]

data <- data_baseline[agecoh_y >= 18]

data[, dur_first_obs_y := dur_first_obs/12]
data[, dur_last_obs_y := dur_last_obs/12]
range(data[, agecoh_y])
range(data[, dur_last_obs_y])

data <- data[dur_first_obs_y >= 0]
data[is.na(marriage), marriage := 0]
data[is.na(separation), separation := 0]

data <- data[!is.na(coh_type)]
data[, coh_type := factor(coh_type, levels = c("prelude", "notready", "conformist",
                                               "reject", "irrelevant"))]
describe(data[, coh_type])

data[, sex := factor(sex_gen, levels = 1:2, labels = c("Men", "Women"))]
data[, east := factor(east, levels = c(1, 0), labels = c("East", "West"))]
data[, educ_cat := factor(educ_cat, levels = c("Low", "Medium", "High"))]
data[, empl := factor(empl)]
data[, prevmar := factor(prevmar, labels = c("No", "Yes"))]
data[, kids2 := factor(ifelse(kids == "child_w_partner", 1, 0),
                       labels = c("no_joint_child", "joint_child"))]
# we have more individuals but numbers are comparable
ftable(data$east, data$coh_type, data$censoring_ind)

# ---- Two time scales analysis same structure as Hiekel et al., 2015

du <- ds <- .5 # half-year
data2d_mar <- prepare_data(
  data = data,
  u = "agecoh_y",
  s_in = "dur_first_obs_y",
  s_out = "dur_last_obs_y",
  events = "marriage",
  covs  = c("sex_gen", "east", "coh_type", "empl", "educ_cat", "kids2", "prevmar", "rel"),
  individual = TRUE,
  ds = ds,
  du = du,
  min_u = 18
)

mod_mar <- fit2ts(data2d_mar,
                  Bbases_spec = list(
                    nseg_u = 11,
                    nseg_s = 12,
                    min_u = 18, max_u = 39,
                    min_s = 0, max_s = 22
                  ),
                  optim_method = "LMMsolver"
)

data2d_sep <- prepare_data(
  data = data,
  u = "agecoh_y",
  s_in = "dur_first_obs_y",
  s_out = "dur_last_obs_y",
  events = "separation",
  covs  = c("sex_gen", "east", "coh_type", "empl", "educ_cat", "kids2", "prevmar", "rel"),
  individual = TRUE,
  ds = ds,
  du = du,
  min_u = 18
)

mod_sep <- fit2ts(data2d_sep,
                  Bbases_spec = list(
                    nseg_u = 11,
                    nseg_s = 12,
                    min_u = 18, max_u = 39,
                    min_s = 0, max_s = 22
                  ),
                  optim_method = "LMMsolver"
)

# let's add the intervals, the B-splines and their specification
# to the saved results
r_Hiekel <- list(
  mod_mar = mod_mar,
  data2d_mar = data2d_mar,
  mod_sep = mod_sep,
  data2d_sep = data2d_sep
)
save(r_Hiekel, file = "R_asHiekel2015.Rda")


# Summaries
summary(mod_mar)
summary(mod_sep)

# ---- Figure S

par(mfrow = c(1,2),
    mar = c(3, 4, 3.3, 1),
    oma = c(2, 0, 0, 1),
    mgp = c(1.5, .5, 0),
    font.main = 1,
    cex.main = 1.2,
    cex.lab = 1.2)
# marriage
plot(mod_mar,
     plot_grid = list(
       c(18, 39, .2),
       c(0, 22, .2)
     ),
     which_plot = "hazard",
     plot_option = list(
       main = "Baseline hazard of marriage",
       xlab = "Attained age",
       ylab = "Duration of the cohabitation"
     )
)

# marriage
plot(mod_sep,
     plot_grid = list(
       c(18, 39, .2),
       c(0, 22, .2)
     ),
     which_plot = "hazard",
     plot_option = list(
       main = "Baseline hazard of separation",
       xlab = "Attained age",
       ylab = "Duration of the cohabitation"
     )
)
