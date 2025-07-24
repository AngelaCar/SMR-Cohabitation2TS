# ---- Figure S2 -----
#---- Kaplan Meier plot of the marginal event (marriage + separation) ----
# Including number of individual at risk in table by group

#---- Load the necessary libraries and functions ----

library(data.table)
library(survival)
library(survminer)

# ---- Load the data ----
load("one_cohabitation18.Rda")

# Indicator for any event
one_cohab[, any_event := ifelse(event_type == "no_event", 0, 1)]

# Variable for 4 groups
one_cohab[sex == 1 & east == 0, group := "MenWest"]
one_cohab[sex == 1 & east == 1, group := "MenEast"]
one_cohab[sex == 2 & east == 0, group := "WomWest"]
one_cohab[sex == 2 & east == 1, group := "WomEast"]
one_cohab[, group := factor(group, 
                            levels = c("WomWest", "WomEast", "MenWest", "MenEast"))]

# ---- KM fit by strata and number of individuals at risk every 5 years ----
km_dur <- survfit(Surv(dur_entry, durcoh_y, any_event) ~ group, 
                  data = one_cohab)

# ---- Plot ----
n1 <- km_dur$strata[1]
n2 <- n1 + km_dur$strata[2]
n3 <- n2 + km_dur$strata[3]
n4 <- n3 + km_dur$strata[4]

par(mfrow = c(1,2),
    mar = c(3, 4, 3.3, .5),
    oma = c(2, 0, 0, 0.5),
    mgp = c(1.5, .5, 0),
    font.main = 1,
    cex.main = 1.2,
    cex.lab = 1.2
)
plot(km_dur$time[1:n1], km_dur$n.risk[1:n1], 
     col = "orchid", pch = 20,
     xlab = "Duration of the cohabitation (years)",
     ylab = "No. of individuals at risk")

points(km_dur$time[n1+1:n2], km_dur$n.risk[n1+1:n2], 
       col = "hotpink", pch = 20)
points(km_dur$time[n2+1:n3], km_dur$n.risk[n2+1:n3], 
       col = "steelblue", pch = 20)
points(km_dur$time[n3+1:n4], km_dur$n.risk[n3+1:n4], 
       col = "turquoise", pch = 20)
legend("topright", pch = 20, col = c("orchid", "hotpink", "steelblue", "turquoise"),
       bty = "n",  c("Women West", "Women East", "Men West", "Men East") )

plot(km_dur$time[1:n1], km_dur$n.event[1:n1], 
     col = "orchid", pch = 20,
     xlab = "Duration of the cohabitation (years)",
     ylab = "No. of events")

points(km_dur$time[n1+1:n2], km_dur$n.event[n1+1:n2], 
       col = "hotpink", pch = 20)
points(km_dur$time[n2+1:n3], km_dur$n.event[n2+1:n3], 
       col = "steelblue", pch = 20)
points(km_dur$time[n3+1:n4], km_dur$n.event[n3+1:n4], 
       col = "turquoise", pch = 20)
abline(v=23, lwd = 0.5, lty = 2, col = "grey50")
legend("topright", pch = 20, col = c("orchid", "hotpink", "steelblue", "turquoise"),
       bty = "n",  c("Women West", "Women East", "Men West", "Men East") )
