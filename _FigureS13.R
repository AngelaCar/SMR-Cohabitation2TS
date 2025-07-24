# ---- Figure S13 ----
# Assessment of proportionality
# ---- Load the necessary libraries and functions ----
library(TwoTimeScales)
library(fields)
library(colorspace)


# ---- Load stratified models ----
load("R_int2d_marriage_me.Rda")
load("R_int2d_marriage_mw.Rda")
load("R_int2d_marriage_we.Rda")
load("R_int2d_marriage_ww.Rda")

# ---- Calculate log-hazard surfaces on a common grid ----
lh_mw <- get_hazard_2d(r_mar_nocov_18p_mw$mod_no_cov,
                       plot_grid = list(
                         c(umin = 18, umax = 47, du = .5),
                         c(smin = 0, smax = 25, ds = .5)
                       ))

lh_ww <- get_hazard_2d(r_mar_nocov_18p_ww$mod_no_cov,
                       plot_grid = list(
                         c(umin = 18, umax = 47, du = .5),
                         c(smin = 0, smax = 25, ds = .5)
                       ))
lh_me <- get_hazard_2d(r_mar_nocov_18p_me$mod_no_cov,
                       plot_grid = list(
                         c(umin = 18, umax = 47, du = .5),
                         c(smin = 0, smax = 25, ds = .5)
                       ))

lh_we <- get_hazard_2d(r_mar_nocov_18p_we$mod_no_cov,
                       plot_grid = list(
                         c(umin = 18, umax = 47, du = .5),
                         c(smin = 0, smax = 25, ds = .5)
                       ))
# common grid
s <- lh_mw$new_plot_grid$ints
u <- lh_mw$new_plot_grid$intu

# ---- Check one: men vs. women, West Germany ----
diff_lh_mw_ww <- lh_mw$loghazard - lh_ww$loghazard
HR_mw_ww <- exp(diff_lh_mw_ww)

# ---- Check two: men vs. women, East Germany ----
diff_lh_me_we <- lh_me$loghazard - lh_we$loghazard
HR_me_we <- exp(diff_lh_me_we)

# ---- Check three: West vs. East, Men ----
diff_lh_mw_me <- lh_mw$loghazard - lh_me$loghazard
HR_mw_me <- exp(diff_lh_mw_me)

# ---- Check four: West vs. East, Women ----
diff_lh_ww_we <- lh_ww$loghazard - lh_we$loghazard
HR_ww_we <- exp(diff_lh_ww_we)

range(HR_mw_ww)
range(HR_me_we)
range(HR_mw_me)
range(HR_ww_we)

comm_breaks <- seq(0.3, 6.4, by = .1)

# HR image plot
col_HR <- diverging_hcl(61, palette = "Blue-Red2", c=70, l=c(60,100))

mat <- matrix(c(1, 2, 3, 4), ncol = 2, nrow = 2, byrow = T)

par(
  mar = c(3, 4, 3.3, .5),
  oma = c(2, 0, 0, 0),
  mgp = c(1.5, .5, 0),
  font.main = 1,
  cex.main = 1.2,
  cex.lab = 1.2
)
layout(mat)

par(mfrow = c(2,2))
image.plot(u,s,HR_mw_ww,
           col = col_HR,
           breaks = comm_breaks,
           main = "Men - Women: West Germany",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(u, s, HR_mw_ww, add=TRUE, 
        col="grey60", labcex=0.8)

image.plot(u,s,HR_me_we,
           col = col_HR,
           breaks = comm_breaks,
           main = "Men - Women: East Germany",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(u, s, HR_me_we, add=TRUE, 
        col="grey60", labcex=0.8)

image.plot(u,s,HR_mw_me,
           col = col_HR,
           breaks = comm_breaks,
           main = "West - East: Men",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(u, s, HR_mw_me, add=TRUE, 
        col="grey60", labcex=0.8)

image.plot(u,s,HR_ww_we,
           col = col_HR,
           breaks = comm_breaks,
           main = "West - East: Women",
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation"
)
contour(u, s, HR_ww_we, add=TRUE, 
        col="grey60", labcex=0.8)
