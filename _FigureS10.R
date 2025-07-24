# ---- Figure S 10  ----
# Hazard of separation -  first cohabitation only ----
# Results of the tts model for the hazard of separation represented in 2d
# in the original plane

# ---- Libraries ----
library(data.table)
library(TwoTimeScales)
library(viridis)

# ---- Load the data ----
load("one_cohabitation18.Rda")

# ---- Restriction to first only ----
first_cohab <- one_cohab[n_cohabitation == 1]

## Separate by gender and east/west
WomWest <- first_cohab[sex == 2 & east == 0]
WomEast <- first_cohab[sex == 2 & east == 1]
MenWest <- first_cohab[sex == 1 & east == 0]
MenEast <- first_cohab[sex == 1 & east == 1]

# identify maximum age at events for the four subgroups
max_age_ww <- WomWest[, max(ageobsend_y)]
max_age_we <- WomEast[, max(ageobsend_y)]
max_age_mw <- MenWest[, max(ageobsend_y)]
max_age_me <- MenEast[, max(ageobsend_y)]

# ---- Load results 2ts models ----
load("R_separation_nocov_firstonly_18p_me.Rda")
load("R_separation_nocov_firstonly_18p_mw.Rda")
load("R_separation_nocov_firstonly_18p_we.Rda")
load("R_separation_nocov_firstonly_18p_ww.Rda")

# ---- Figure ----
pal_2ts <- function(n, alpha) {
  rev(plasma(n, alpha))
}
mat <- matrix(c(1, 2, 3, 4), ncol = 2, nrow = 2, byrow = T)
breaks <- seq(0, 0.35, length = 51)

par(
  mar = c(3, 4, 3.3, .5),
  oma = c(2, 0, 0, 0.5),
  mgp = c(1.5, .5, 0),
  font.main = 1,
  cex.main = 1.2,
  cex.lab = 1.2
)
layout(mat)

# Men East
plot(x = r_sep_nocov_frst_me$mod_no_cov,
     plot_grid = list(
       c(18, 46, .2),
       c(0, 26, .2)
     ),
     which_plot = "hazard",
     plot_option = list(
       original = T,
       rectangular_grid = F,
       tmax = max_age_me,
       col_palette = pal_2ts,
       breaks = breaks,
       main = "Men East",
       xlab = "Attained age",
       ylab = "Duration of the cohabitation"
     )
)
polygon(
  x = c(max_age_me - .2, max_age_me - .2, max_age_me + ((max_age_me / 100) * 2), max_age_me + ((max_age_me / 100) * 2)),
  y = c(0, 26.5, 26.5, 0),
  col = "white", border = NA
)

plot(x = r_sep_nocov_frst_mw$mod_no_cov,
     plot_grid = list(
       c(18, 47, .2),
       c(0, 25, .2)
     ),
     which_plot = "hazard",
     plot_option = list(
       original = T,
       rectangular_grid = F,
       tmax = max_age_mw,
       col_palette = pal_2ts,
       breaks = breaks,
       main = "Men West",
       xlab = "Attained age",
       ylab = "Duration of the cohabitation"
     )
)
polygon(
  x = c(max_age_mw - .2, max_age_mw - .2, max_age_mw + ((max_age_mw / 100) * 2),
        max_age_mw + ((max_age_mw / 100) * 2)),
  y = c(0, 25.5, 25.5, 0),
  col = "white", border = NA
)

plot(x = r_sep_nocov_frst_we$mod_no_cov,
     plot_grid = list(
       c(18, 43, .2),
       c(0, 25, .2)
     ),
     which_plot = "hazard",
     plot_option = list(
       original = T,
       rectangular_grid = F,
       tmax = max_age_we,
       col_palette = pal_2ts,
       breaks = breaks,
       main = "Women East",
       xlab = "Attained age",
       ylab = "Duration of the cohabitation"
     )
)
polygon(
  x = c(max_age_we - .2, max_age_we - .2, max_age_we + ((max_age_we / 100) * 2), 
        max_age_we + ((max_age_we / 100) * 2)),
  y = c(0, 25.5, 25.5, 0),
  col = "white", border = NA
)

plot(x = r_sep_nocov_frst_ww$mod_no_cov,
     plot_grid = list(
       c(18, 44, .2),
       c(0, 27, .2)
     ),
     which_plot = "hazard",
     plot_option = list(
       original = T,
       rectangular_grid = F,
       tmax = max_age_ww,
       col_palette = pal_2ts,
       breaks = breaks,
       main = "Women West",
       xlab = "Attained age",
       ylab = "Duration of the cohabitation"
     )
)
polygon(
  x = c(max_age_ww - .2, max_age_ww - .2, max_age_ww + ((max_age_ww / 100) * 2),
        max_age_ww + ((max_age_ww / 100) * 2)),
  y = c(0, 27.5, 27.5, 0),
  col = "white", border = NA
)


