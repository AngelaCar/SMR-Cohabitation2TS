#############################################################
## Analysis of transitions out of cohabitation (pairfam)   ##
## Data preparation                                        ##
## Angela Carollo                                          ##
#############################################################

# ---- libraries and stata dataset ----
library(readstata13) # for reading the Stata files
library(Hmisc)       # for function `describe`
library(data.table)  # for easy manipulation

# ----- data cleaning ----- 

# load the data - here please change the path to your data folder!!!
biopart <- read.dta13("K:/Pairfam/Release11-0/Data/Stata/biopart.dta",
                      convert.factors = FALSE,
                      replace.strl = TRUE)


# labels of the variables as in stata
labels <- get.label.tables(biopart)

# include only cohabitation
biopart_coh <- subset(biopart, cohbeg > 0) 
# positive numbers indicate dates of beginning of 
# the cohabitation

# exclude same-sex couples
biopart_het_coh <- subset(biopart_coh, homosex_p == 0)
# observations excluded: 172 (not individuals but relationships)
# nrow(biopart_coh) - nrow(biopart_het_coh)

# Number of cohabitation for individual
biopart_het_coh <- as.data.table(biopart_het_coh)
biopart_het_coh[, n_cohabitation := 1:.N, by = id]
biopart_het_coh[, id_rel := 1:.N]
biopart_het_coh <- as.data.frame(biopart_het_coh)

# Which columns have the interview dates?
intdatw <- c("intdatw1", "intdatw2", "intdatw3", "intdatw4", "intdatw5",
             "intdatw6", "intdatw7", "intdatw8", "intdatw9", "intdatw10",
             "intdatw11")
wcol_intd <- which(names(biopart_het_coh) %in% intdatw)

# Date of last interview per individual
for(i in 1:nrow(biopart_het_coh)){
  biopart_het_coh$lastinterview[i] <- max(biopart_het_coh[i, wcol_intd])
  
}
# describe(biopart_het_coh$lastinterview)
# will be useful to impute end of observation period for
# censored observations

# Find the wave just after the cohabitation started.
# this is to link covariates from the wave's specific dataset to the time after cohabitation

# wave after cohabitation
for (i in 1:nrow(biopart_het_coh)) {
  biopart_het_coh$wave_Acoh[i] <- ifelse(length(which(biopart_het_coh[i, wcol_intd] >= biopart_het_coh$cohbeg[i])) != 0,
    biopart_het_coh[i, wcol_intd[min(which(biopart_het_coh[i, wcol_intd] >= biopart_het_coh$cohbeg[i]))]],
    NA
  )
}
# describe(biopart_het_coh$wave_Acoh)

# wave before cohabitation
for (i in 1:nrow(biopart_het_coh)) {
  biopart_het_coh$wave_Bcoh[i] <- ifelse(length(which((biopart_het_coh[i, wcol_intd] < biopart_het_coh$cohbeg[i] ) & (biopart_het_coh[i, wcol_intd] > 0))) != 0,
    biopart_het_coh[i, wcol_intd[max(which((biopart_het_coh[i, wcol_intd] < biopart_het_coh$cohbeg[i]) & (biopart_het_coh[i, wcol_intd] > 0)))]],
    NA
  )
}
# describe(biopart_het_coh$wave_Bcoh)

# there was one case where the beginning of the cohabitation started one month 
# after the interview date and then the individual dropped out of the sample,
# therefore we drop this observation from the sample too.
biopart_het_coh <- subset(biopart_het_coh, !is.na(wave_Acoh))

# Coding cohabitation status by combining info from cohend and marbeg and marend
table(biopart_het_coh$marbeg == -7)
biopart_het_coh <- subset(biopart_het_coh, marbeg != -7)
biopart_het_coh <- as.data.table(biopart_het_coh)
table(biopart_het_coh$cohend == -66)
biopart_het_coh <- subset(biopart_het_coh, cohend != -66)

biopart_het_coh[, status := ifelse(cohend == -99 & marbeg == -3, "coh_ongoing", NA)]
biopart_het_coh[, status := ifelse(cohend == -99 & marbeg > 0 & marend == -99 | marend == -7, "mar_ongoing", status)]
biopart_het_coh[, status := ifelse(cohend == -99 & marbeg > 0 & marend > 0, "div_cohabiting", status)]
biopart_het_coh[, status := ifelse(cohend > 0 & marbeg > 0 & marend == -99, "mar_separated", status)]
biopart_het_coh[, status := ifelse(cohend > 0 & marbeg > 0 & marend > 0, "mar_divorced", status)]
biopart_het_coh[, status := ifelse(cohend > 0 & marbeg == -3, "coh_dissolved", status)]
biopart_het_coh[is.na(status), status := "unclear"]

describe(biopart_het_coh$status)

# remove unclear cases
biopart_het_coh <- biopart_het_coh[ status != "unclear"]

# Ending of cohabitation according to status
biopart_het_coh[status == "coh_dissolved", cohend_new := cohend]
biopart_het_coh[status %in% c("mar_ongoing", "mar_divorced", "mar_separated",
                              "div_cohabiting") , cohend_new := marbeg]
biopart_het_coh[status == "coh_ongoing", cohend_new := -99]

biopart_het_coh <- as.data.frame(biopart_het_coh)

# First interview 
for(i in 1:nrow(biopart_het_coh)){
  biopart_het_coh$firstinterview[i] <- biopart_het_coh[i, wcol_intd[min(which(biopart_het_coh[i, wcol_intd] > 0))]] 
}

# Select only cohabitation which did not end before the first interview

for(i in 1:nrow(biopart_het_coh)){
  biopart_het_coh$ended_before_first[i] <- ifelse( 
    biopart_het_coh$firstinterview[i] >= (biopart_het_coh$cohend_new[i]) & (biopart_het_coh$cohend_new[i] > 0),
    1, 0)
}
describe(biopart_het_coh$ended_before_first)

with(subset(biopart_het_coh, ended_before_first == 1), table(status))
biopart_het_coh2 <- biopart_het_coh[biopart_het_coh$ended_before_first == 0,]

# keep only a subset of covariates
covariates <- c("id", "demodiff", "cohort", "sex", "dob", 
                "partindex", "sexp", "dobp", "dodp", 
                "cohbeg", "cohend",
                "marbeg", "marend", "marcer", 
                "intdatw1", "intdatw2", "intdatw3", "intdatw4", "intdatw5", "intdatw6",
                "intdatw7", "intdatw8", "intdatw9", "intdatw10", "intdatw11", "intdatw12",
                "firstinterview",
                "lastinterview", "wave_Acoh", "wave_Bcoh", "n_cohabitation", "cohend_new", "status")
biopart_het_coh2 <- biopart_het_coh2[, names(biopart_het_coh2) %in% covariates]

biopart_het_coh2 <- as.data.table(biopart_het_coh2)

# has the cohabitation started after a marriage or with the marriage?
biopart_het_coh2[, cohabaftermarr := ifelse((marbeg >= 0 & cohbeg >= marbeg), 1, 0)]
# table(biopart_het_coh2$cohabaftermarr)

# remove all cohabitation which started after/with a marriage
biopart_het_coh2 <- biopart_het_coh2[ cohabaftermarr == 0]

# negative values indicate a censored observation
# with(subset(biopart_het_coh2, cohend_new < 0), table(cohend_new))
# -99 is a censored observation, which includes both individuals who dropped out of the sample
# and individuals who did not marry yet.
# For those, the date out of cohabitation will be the date of last interview
biopart_het_coh2[, cohend_new := ifelse((cohend_new == -99 & lastinterview >= cohbeg), 
                                           lastinterview, cohend_new)]
describe(biopart_het_coh2$status)

# Event type
biopart_het_coh2[, event_type := ifelse(cohend_new == marbeg, "marriage", NA)]
biopart_het_coh2[, event_type := ifelse(status == "coh_dissolved", "separation", event_type)]
biopart_het_coh2[, event_type := ifelse(status == "coh_ongoing", "no_event", event_type)]
# describe(biopart_het_coh2$event_type)

# Time of entry in our sample (maximum between time of cohabitation
# and first interview)
biopart_het_coh2[, entry_time := pmax(firstinterview, cohbeg)]
#View(biopart_het_coh2[, c("id", "wave_coh", "cohbeg", "entry_time")])

# Marriage begins before first entry time?
biopart_het_coh2[, marbeforeentry := ifelse(event_type == "marriage" & (marbeg < entry_time), 1, 0)]
# table(biopart_het_coh2$marbeforeentry)
biopart_het_coh2[, sepbeforeentry := ifelse(event_type == "separation" & cohend_new < entry_time, 1, 0)]
# table(biopart_het_coh2$sepbeforeentry)

# age variables (in months and years)
biopart_het_coh2[, agecoh := (cohbeg - dob)]
biopart_het_coh2[, ageobsend := (cohend_new - dob)]
biopart_het_coh2[, durcoh := (cohend_new - cohbeg)]

biopart_het_coh2[, agecoh_y := (cohbeg - dob)/12]
biopart_het_coh2[, ageobsend_y := (cohend_new - dob)/12]
biopart_het_coh2[, durcoh_y := (cohend_new - cohbeg)/12]
# summary(biopart_het_coh2$agecoh)
# summary(biopart_het_coh2$ageobsend)
# summary(biopart_het_coh2$durcoh)
# table(biopart_het_coh2$agecoh_y < 18)

# age restrictions on start of the cohabitation 
biopart_het_coh2 <- biopart_het_coh2[agecoh_y >= 18]
# restrictions on marriages below 18 years
biopart_het_coh2 <- biopart_het_coh2[!((ageobsend_y < 18) & (event_type == "marriage"))]

# ---- Covariates ----
# Include variables which count the cohabitation
# information about living in east or west germany
# with interview date
all_waves <- read.dta13("anchor_all_waves.dta",
                        convert.factors = FALSE,
                        replace.strl = TRUE)
all_waves <- as.data.table(all_waves)

indiv <- biopart_het_coh2[!duplicated(biopart_het_coh2$id),]$id
all_waves_ind <- all_waves[id %in% indiv]

# merge the data by interview date and wave_coh
biopart_het_coh3 <- merge(biopart_het_coh2, all_waves_ind, all.x = TRUE,  
                         by.x = c("id", "wave_Acoh"), by.y = c("id", "dateint"))

# as.data.table
biopart_het_coh3 <- as.data.table(biopart_het_coh3)

# order by id and beginning of cohabitation
setorder(biopart_het_coh3, id, cohbeg)

# Education
education <- all_waves_ind[, c("id","ageint_y", "isced", "isced2")]
#describe(biopart_het_coh3$isced)
#describe(biopart_het_coh3$isced2)

biopart_het_coh3[, educ_cat := ifelse(isced2 == -7 | isced2 == 0, "Unknown", NA)]
biopart_het_coh3[, educ_cat := ifelse(isced2 == 1 | isced2 == 2 | isced2 == 3,
                                 "Low", educ_cat)]
biopart_het_coh3[, educ_cat := ifelse(isced2 == 4 | isced2 == 5 | isced2 == 6,
                                 "Middle", educ_cat)]
biopart_het_coh3[, educ_cat := ifelse(isced2 == 7 | isced2 == 8 ,
                                 "High", educ_cat)]
describe(biopart_het_coh3$educ_cat)

# remove Unknown education
biopart_het_coh4 <- biopart_het_coh3[!educ_cat == "Unknown"]

table(biopart_het_coh4$east)
# one missing values for east
biopart_het_coh4 <- biopart_het_coh4[! east == -7]
table(biopart_het_coh4$cohort)
# only one individual from cohort 4
biopart_het_coh4 <- biopart_het_coh4[! cohort == 4]
#

# rename the data
all_cohab_waves <- biopart_het_coh4

# SAVE THE DATA --- these are ALL cohabitations, so no restriction to first 
# cohabitation yet.
# save(all_cohab_waves, file = "all_cohabit_w_cov_1-11.Rda")

#---- First cohabitation data prep ----
# load("all_cohabit_w_cov_1-11.Rda")
## One cohabitation per individual only
setorder(all_cohab_waves, id, cohbeg)
# table(all_cohab_waves$n_cohabitation)
describe(all_cohab_waves$id)
one_cohab <- all_cohab_waves[!duplicated(id)]

# Duration of cohabitation at entry time (entry_time - cohbeg)/12 (in years)
one_cohab[, dur_entry := ((entry_time - cohbeg)/12)]
one_cohab[, age_entry := agecoh_y + dur_entry]
# Is this the first cohabitation ever or not?
one_cohab[, first_ever := ifelse(n_cohabitation == 1, "yes", "no")]

# finally, censor events that happen after more than 23 years of cohabitation
# there are only 3
View(one_cohab[durcoh_y >= 23 & event_type != "no_event"])
one_cohab[durcoh_y >= 23 & event_type != "no_event", event_type := "no_event"]

save(one_cohab, file = "one_cohabitation18.Rda")


## For sensitivity analysis, save one version of the file where age at cohabitation
## is >= 15
# biopart_het_coh2_15p <- biopart_het_coh2[agecoh_y >= 15]
# # restrictions on marriages below 18 years
# biopart_het_coh2_15p <- biopart_het_coh2_15p[!((ageobsend_y < 18) & (event_type == "marriage"))]

# merge the data by interview date and wave_coh
# biopart_het_coh3_15p <- merge(biopart_het_coh2_15p, all_waves_ind, all.x = TRUE,  
#                           by.x = c("id", "wave_Acoh"), by.y = c("id", "dateint"))
# 
# # as.data.table
# biopart_het_coh3_15p <- as.data.table(biopart_het_coh3_15p)
# 
# # order by id and beginning of cohabitation
# setorder(biopart_het_coh3_15p, id, cohbeg)
# 
# # Education
# education <- all_waves_ind[, c("id","ageint_y", "isced", "isced2")]
# #describe(biopart_het_coh3$isced)
# #describe(biopart_het_coh3$isced2)
# 
# biopart_het_coh3_15p[, educ_cat := ifelse(isced2 == -7 | isced2 == 0, "Unknown", NA)]
# biopart_het_coh3_15p[, educ_cat := ifelse(isced2 == 1 | isced2 == 2 | isced2 == 3,
#                                       "Low", educ_cat)]
# biopart_het_coh3_15p[, educ_cat := ifelse(isced2 == 4 | isced2 == 5 | isced2 == 6,
#                                       "Middle", educ_cat)]
# biopart_het_coh3_15p[, educ_cat := ifelse(isced2 == 7 | isced2 == 8 ,
#                                       "High", educ_cat)]
# describe(biopart_het_coh3_15p$educ_cat)
# 
# # remove Unknown education
# biopart_het_coh4_15p <- biopart_het_coh3_15p[!educ_cat == "Unknown"]
# 
# table(biopart_het_coh4_15p$east)
# # one missing values for east
# biopart_het_coh4_15p <- biopart_het_coh4_15p[! east == -7]
# table(biopart_het_coh4_15p$cohort)
# # only one individual from cohort 4
# biopart_het_coh4_15p <- biopart_het_coh4_15p[! cohort == 4]
# #
# 
# # rename the data
# all_cohab_waves_15p <- biopart_het_coh4_15p
# 
# ## One cohabitation per individual only
# setorder(all_cohab_waves_15p, id, cohbeg)
# # table(all_cohab_waves$n_cohabitation)
# describe(all_cohab_waves_15p$id)
# one_cohab_15p <- all_cohab_waves_15p[!duplicated(id)]
# 
# # Duration of cohabitation at entry time (entry_time - cohbeg)/12 (in years)
# one_cohab_15p[, dur_entry := ((entry_time - cohbeg)/12)]
# one_cohab_15p[, age_entry := agecoh_y + dur_entry]
# # Is this the first cohabitation ever or not?
# one_cohab_15p[, first_ever := ifelse(n_cohabitation == 1, "yes", "no")]
# 
# # finally, censor events that happen after more than 23 years of cohabitation
# # there are only 3
# View(one_cohab_15p[durcoh_y >= 23 & event_type != "no_event"])
# one_cohab_15p[durcoh_y >= 23 & event_type != "no_event", event_type := "no_event"]
# 
# save(one_cohab_15p, file = "one_cohabitation15.Rda")

