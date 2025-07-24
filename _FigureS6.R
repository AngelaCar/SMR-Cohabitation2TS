# ---- Figure S6 -----
# Hazard of marriage over duration of the cohabitation
# selected values of the age at entry into cohabitation
# comparison of 3 models for Women East and Women West

#
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# ------------------------------------------------
rm(list=ls())
# ------------------------------------------------

library(TwoTimeScales)
library(JOPS)

# ---- Results ----
load("R_marriage_vcm_we.Rda")
load("R_marriage_vcm_ww.Rda")
load("R_int2d_marriage_we.Rda")
load("R_int2d_marriage_ww.Rda")
load("R_pgam_marriage_we.Rda")
load("R_pgam_marriage_ww.Rda")

# ----- Estimate hazard with denser grid -----
# First: make denser grid
# Second: re-evaluate B-splines
# Third: compute hazards

# Men East
bins_we <- r_pgam_marr_we$bins 
we_bins_u_new <- seq(20, 45, by=5)
we_bins_s_new <- seq(min(bins_we$bins_s), max(bins_we$bins_s), by = .2)

# For Pgam model: these new intervals should be expanded to a grid
we_new_grid <- expand.grid(we_bins_u_new, we_bins_s_new)
names(we_new_grid) <- c("bins_u_new", "bins_s_new")

we_pgam_Bu_new <- JOPS::bbase(we_new_grid$bins_u_new, 
                              xl = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                              xr = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                              nseg = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

we_pgam_Bs_new <- JOPS::bbase(we_new_grid$bins_s_new, 
                              xl = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                              xr = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                              nseg = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

we_pgam_B_new <- cbind(we_pgam_Bu_new, we_pgam_Bs_new)

# For VCM and int2D we only need the marginal B-splines
we_Bu_new <- JOPS::bbase(we_bins_u_new, 
                         xl = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                         xr = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                         nseg = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

we_Bs_new <- JOPS::bbase(we_bins_s_new, 
                         xl = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                         xr = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                         nseg = attr(r_mar_nocov_18p_we$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

# Men West
bins_ww <- r_pgam_marr_ww$bins 
ww_bins_u_new <- seq(20, 45, by= 5)
ww_bins_s_new <- seq(min(bins_ww$bins_s), max(bins_ww$bins_s), by = .2)

# For Pgam model: these new intervals should be expanded to a grid
ww_new_grid <- expand.grid(ww_bins_u_new, ww_bins_s_new)
names(ww_new_grid) <- c("bins_u_new", "bins_s_new")

ww_pgam_Bu_new <- JOPS::bbase(ww_new_grid$bins_u_new, 
                              xl = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                              xr = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                              nseg = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

ww_pgam_Bs_new <- JOPS::bbase(ww_new_grid$bins_s_new, 
                              xl = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                              xr = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                              nseg = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

ww_pgam_B_new <- cbind(ww_pgam_Bu_new, ww_pgam_Bs_new)

# For VCM and int2D we only need the marginal B-splines
ww_Bu_new <- JOPS::bbase(ww_bins_u_new, 
                         xl = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                         xr = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                         nseg = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

ww_Bs_new <- JOPS::bbase(ww_bins_s_new, 
                         xl = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                         xr = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                         nseg = attr(r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

# Calculate hazards
#  Pgam
Eta_pgam_we <- we_pgam_B_new %*% r_pgam_marr_we$mod_pgam_marr_we$alpha
Haz_pgam_we <- exp(Eta_pgam_we)
Haz_pgam_we <- matrix(Haz_pgam_we, length(we_bins_u_new), length(we_bins_s_new))

Eta_pgam_ww <- ww_pgam_B_new %*% r_pgam_marr_ww$mod_pgam_marr_ww$alpha
Haz_pgam_ww <- exp(Eta_pgam_ww)
Haz_pgam_ww <- matrix(Haz_pgam_ww, length(ww_bins_u_new), length(ww_bins_s_new))

# VCM
eta_we <- we_Bs_new %*% r_mar_vcm_we$vcm_mar_we$Alpha[, 1]
vc_we  <- we_Bs_new %*% r_mar_vcm_we$vcm_mar_we$Alpha[, 2]
eta_ww <- ww_Bs_new %*% r_mar_vcm_ww$vcm_mar_ww$Alpha[, 1]
vc_ww  <- ww_Bs_new %*% r_mar_vcm_ww$vcm_mar_ww$Alpha[, 2]

HazVCM_we <- exp( matrix(rep(eta_we, length(we_bins_u_new)), 
                         byrow = T, nrow = length(we_bins_u_new), 
                         ncol = length(we_bins_s_new)) + 
                    outer(we_bins_u_new, as.vector(vc_we))  )
HazVCM_ww <- exp( matrix(rep(eta_ww, length(ww_bins_u_new)), 
                         byrow = T, nrow = length(ww_bins_u_new), 
                         ncol = length(ww_bins_s_new)) + 
                    outer(ww_bins_u_new, as.vector(vc_ww))  )

#  2D interaction
Haz2D_we <- exp(we_Bu_new %*% r_mar_nocov_18p_we$mod_no_cov$optimal_model$Alpha %*% t(we_Bs_new))
Haz2D_ww <- exp(ww_Bu_new %*% r_mar_nocov_18p_ww$mod_no_cov$optimal_model$Alpha %*% t(ww_Bs_new))

# ---- Figure comparing cross-sections of the hazard from the three models ----
legrange <- range(Haz_pgam_we, HazVCM_we, Haz2D_we,
                  Haz_pgam_ww, HazVCM_ww, Haz2D_ww)

mat <-  matrix(1:6, ncol = 3, nrow = 2, byrow = T)
breaks <- seq(0, legrange[2]+.01, length = 100)

par(mar= c(3,4,4,.5),
    oma = c(2, 0, 2,0),
    mgp = c(1.5,.5,0),
    cex.main = 1.5,
    cex.lab = 1.2,
    font.main = 1)
layout(mat)

grey7 <- grey.colors(n = 7)

matplot(x = ww_bins_s_new, y = t(Haz_pgam_ww), 
        type = "l", col = grey.colors(n = 7), #150), 
        lty = 1, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")
matplot(x = ww_bins_s_new, y = t(HazVCM_ww), 
        type = "l", col = grey.colors(n = 7), # = 150), 
        lty = 1, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")
matplot(x = ww_bins_s_new, y = t(Haz2D_ww), 
        type = "l", col = grey.colors(n = 7), #150), 
        lty = 1, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")


matplot(x = we_bins_s_new, y = t(Haz_pgam_we), 
        type = "l", 
        col = grey.colors(n = 7), # 150), 
        lty = 1, 
        ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")

# ---- legend goes here (bootom row, leftmost)
legend(x=0, y=0.25, ncol =2, bty="n", cex = 1.0,
       legend = paste(" u =", seq(20, 45, by=5)), lty=1, col=grey7)


matplot(x = we_bins_s_new, y = t(HazVCM_we), 
        type = "l", col = grey.colors(n =7), # = 150), 
        lty = 1, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")
matplot(x = we_bins_s_new, y = t(Haz2D_we), 
        type = "l", col = grey.colors(n = 7), #150), 
        lty = 1, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")


# Add common title for the first row
mtext("Women West Germany", 
      side = 3, line = -3.4, at = 0.52,
      outer = TRUE, cex = 1.25)
# Common title for second row
mtext("Women East Germany", 
      side = 3, line = -23.9, at = 0.52,
      outer = TRUE, cex = 1.25)

# Column titles
mtext("Model (A):\n log-additive", 
      side = 3, line = -0.9, at = 0.18, 
      outer = TRUE, cex = 1.0)
mtext("Model (B):\n varying-coefficients", 
      side = 3, line = -0.9, at = 0.53, 
      outer = TRUE, cex = 1.0)
mtext("Model (C):\n 2D interaction", 
      side = 3, line = -0.9, at = 0.86, 
      outer = TRUE, cex = 1.0)
