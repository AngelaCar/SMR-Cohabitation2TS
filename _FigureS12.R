# ---- Figure S12 ----
# Plot of SEs for the marriage model and comparison with sample first cohabitation only 

# ---- Libraries ----
library(data.table)
library(TwoTimeScales)
library(colorspace)
library(fields)

# ---- Data ----
load("one_cohabitation18.Rda")

## Separate by gender and east/west
WomWest <- one_cohab[sex == 2 & east == 0]
WomEast <- one_cohab[sex == 2 & east == 1]
MenWest <- one_cohab[sex == 1 & east == 0]
MenEast <- one_cohab[sex == 1 & east == 1]

# identify maximum age at events for the four subgroups
max_age_ww <- WomWest[, max(ageobsend_y)]
max_age_we <- WomEast[, max(ageobsend_y)]
max_age_mw <- MenWest[, max(ageobsend_y)]
max_age_me <- MenEast[, max(ageobsend_y)]

# ---- Load results models age 18+ sample ----
load("R_int2d_marriage_me.Rda")
load("R_int2d_marriage_mw.Rda")
load("R_int2d_marriage_we.Rda")
load("R_int2d_marriage_ww.Rda")

# ---- Load results models age 18+, first cohabitation sample ----
load("R_marriage_nocov_firstonly_18p_me.Rda")
load("R_marriage_nocov_firstonly_18p_mw.Rda")
load("R_marriage_nocov_firstonly_18p_we.Rda")
load("R_marriage_nocov_firstonly_18p_ww.Rda")

# ---- Save SEs in matrices ----
p_SE_full_me <- get_hazard_2d(r_mar_nocov_18p_me$mod_no_cov,
                              plot_grid = list(
                                c(18, 47, .2),
                                c(0, 26, .2)
                              ))

p_SE_full_mw <- get_hazard_2d(r_mar_nocov_18p_mw$mod_no_cov,
                              plot_grid = list(
                                c(18, 47, .2),
                                c(0, 25, .2)
                              ))

p_SE_full_we <- get_hazard_2d(r_mar_nocov_18p_we$mod_no_cov,
                              plot_grid = list(
                                c(18, 47, .2),
                                c(0, 25, .2)
                              ))

p_SE_full_ww <- get_hazard_2d(r_mar_nocov_18p_ww$mod_no_cov,
                              plot_grid = list(
                                c(18, 47, .2),
                                c(0, 27, .2)
                              ))

p_SE_frst_me <- get_hazard_2d(r_mar_nocov_frst_me$mod_no_cov,
                              plot_grid = list(
                                c(18, 46, .2),
                                c(0, 26, .2)
                              ))

p_SE_frst_mw <- get_hazard_2d(r_mar_nocov_frst_mw$mod_no_cov,
                              plot_grid = list(
                                c(18, 47, .2),
                                c(0, 25, .2)
                              ))

p_SE_frst_we <- get_hazard_2d(r_mar_nocov_frst_we$mod_no_cov,
                              plot_grid = list(
                                c(18, 43, .2),
                                c(0, 25, .2)
                              ))

p_SE_frst_ww <- get_hazard_2d(r_mar_nocov_frst_ww$mod_no_cov,
                              plot_grid = list(
                                c(18, 44, .2),
                                c(0, 27, .2)
                              ))

# ---- Figure ----
pal_SE <- rev(sequential_hcl(n = 50, "Red-Purple"))
mat <- matrix(1:8, ncol = 2, nrow = 4, byrow = T)

range(p_SE_full_me$SE_hazard, p_SE_full_mw$SE_hazard, 
      p_SE_full_we$SE_hazard, p_SE_full_ww$SE_hazard,
      p_SE_frst_me$SE_hazard, p_SE_frst_mw$SE_hazard,
      p_SE_frst_we$SE_hazard, p_SE_frst_ww$SE_hazard)
breaks <- seq(0, 0.6, length = 51)

par(
  mar = c(3, 4, 3.3, .5),
  oma = c(2, 0, 0, 0),
  mgp = c(1.5, .5, 0),
  cex.lab = 1.2,
  cex.main = 1.2,
  font.main = 1
)
layout(mat)

image.plot(p_SE_full_me$new_plot_grid$intu, p_SE_full_me$new_plot_grid$ints,
           p_SE_full_me$SE_hazard,
           col = pal_SE,
           breaks = breaks,
           main = "Men East - all cohabitations",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(p_SE_full_me$new_plot_grid$intu, p_SE_full_me$new_plot_grid$ints,
        p_SE_full_me$SE_hazard,
        col = "grey",
        nlevels = 10,
        add = T
)
image.plot(p_SE_frst_me$new_plot_grid$intu, p_SE_frst_me$new_plot_grid$ints,
           p_SE_frst_me$SE_hazard,
           col = pal_SE,
           breaks = breaks,
           main = "Men East - first cohabitations only",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(p_SE_frst_me$new_plot_grid$intu, p_SE_frst_me$new_plot_grid$ints,
        p_SE_frst_me$SE_hazard,
        col = "grey",
        nlevels = 10,
        add = T
)

image.plot(p_SE_full_mw$new_plot_grid$intu, p_SE_full_mw$new_plot_grid$ints,
           p_SE_full_mw$SE_hazard,
           col = pal_SE,
           breaks = breaks,
           main = "Men West - all cohabitations",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(p_SE_full_mw$new_plot_grid$intu, p_SE_full_mw$new_plot_grid$ints,
        p_SE_full_mw$SE_hazard,
        col = "grey",
        nlevels = 10,
        add = T
)
image.plot(p_SE_frst_mw$new_plot_grid$intu, p_SE_frst_mw$new_plot_grid$ints,
           p_SE_frst_mw$SE_hazard,
           col = pal_SE,
           breaks = breaks,
           main = "Men West - first cohabitations only",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(p_SE_frst_mw$new_plot_grid$intu, p_SE_frst_mw$new_plot_grid$ints,
        p_SE_frst_mw$SE_hazard,
        col = "grey",
        nlevels = 10,
        add = T
)

image.plot(p_SE_full_we$new_plot_grid$intu, p_SE_full_we$new_plot_grid$ints,
           p_SE_full_we$SE_hazard,
           col = pal_SE,
           breaks = breaks,
           main = "Women East - all cohabitations",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(p_SE_full_we$new_plot_grid$intu, p_SE_full_we$new_plot_grid$ints,
        p_SE_full_we$SE_hazard,
        col = "grey",
        nlevels = 10,
        add = T
)
image.plot(p_SE_frst_we$new_plot_grid$intu, p_SE_frst_we$new_plot_grid$ints,
           p_SE_frst_we$SE_hazard,
           col = pal_SE,
           breaks = breaks,
           main = "Women East - first cohabitations only",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(p_SE_frst_we$new_plot_grid$intu, p_SE_frst_we$new_plot_grid$ints,
        p_SE_frst_we$SE_hazard,
        col = "grey",
        nlevels = 10,
        add = T
)

image.plot(p_SE_full_ww$new_plot_grid$intu, p_SE_full_ww$new_plot_grid$ints,
           p_SE_full_ww$SE_hazard,
           col = pal_SE,
           breaks = breaks,
           main = "Women West - all cohabitations",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(p_SE_full_ww$new_plot_grid$intu, p_SE_full_ww$new_plot_grid$ints,
        p_SE_full_ww$SE_hazard,
        col = "grey",
        nlevels = 10,
        add = T
)
image.plot(p_SE_frst_ww$new_plot_grid$intu, p_SE_frst_ww$new_plot_grid$ints,
           p_SE_frst_ww$SE_hazard,
           col = pal_SE,
           breaks = breaks,
           main = "Women West - first cohabitations only",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(p_SE_frst_ww$new_plot_grid$intu, p_SE_frst_ww$new_plot_grid$ints,
        p_SE_frst_ww$SE_hazard,
        col = "grey",
        nlevels = 10,
        add = T
)

