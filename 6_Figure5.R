# ---- Figure 5 -----
# Hazard of marriage over duration of the cohabitation
# selected values of the age at entry into cohabitation
# comparison of 3 models for Men East and Men West

#
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# ------------------------------------------------
rm(list=ls())
# ------------------------------------------------


library(fields)
library(viridis)
library(TwoTimeScales)
library(RColorBrewer)
library(JOPS)

# ---- Results ----
load("R_marriage_vcm_me.Rda")
load("R_marriage_vcm_mw.Rda")
load("R_int2d_marriage_me.Rda")
load("R_int2d_marriage_mw.Rda")
load("R_pgam_marriage_me.Rda")
load("R_pgam_marriage_mw.Rda")

# ----- Estimate hazard with denser grid -----
# First: make denser grid
# Second: re-evaluate B-splines
# Third: compute hazards

# Men East
bins_me <- r_mar_pgam_me$bins 
me_bins_u_new <- seq(20, 45, by = 5) #min(bins_me$bins_u), max(bins_me$bins_u), by = .2)
me_bins_s_new <- seq(min(bins_me$bins_s), max(bins_me$bins_s), by = .2)

# For Pgam model: these new intervals should be expanded to a grid
me_new_grid <- expand.grid(me_bins_u_new, me_bins_s_new)
names(me_new_grid) <- c("bins_u_new", "bins_s_new")

me_pgam_Bu_new <- JOPS::bbase(me_new_grid$bins_u_new, 
                              xl = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                              xr = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                              nseg = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

me_pgam_Bs_new <- JOPS::bbase(me_new_grid$bins_s_new, 
                              xl = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                              xr = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                              nseg = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

me_pgam_B_new <- cbind(me_pgam_Bu_new, me_pgam_Bs_new)

# For VCM and int2D we only need the marginal B-splines
me_Bu_new <- JOPS::bbase(me_bins_u_new, 
                         xl = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                         xr = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                         nseg = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

me_Bs_new <- JOPS::bbase(me_bins_s_new, 
                         xl = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                         xr = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                         nseg = attr(r_mar_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

# Men West
bins_mw <- r_mar_pgam_mw$bins 
mw_bins_u_new <- seq(20, 45, by =  5)
mw_bins_s_new <- seq(min(bins_mw$bins_s), max(bins_mw$bins_s), by = .2)

# For Pgam model: these new intervals should be expanded to a grid
mw_new_grid <- expand.grid(mw_bins_u_new, mw_bins_s_new)
names(mw_new_grid) <- c("bins_u_new", "bins_s_new")

mw_pgam_Bu_new <- JOPS::bbase(mw_new_grid$bins_u_new, 
                              xl = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                              xr = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                              nseg = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

mw_pgam_Bs_new <- JOPS::bbase(mw_new_grid$bins_s_new, 
                              xl = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                              xr = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                              nseg = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

mw_pgam_B_new <- cbind(mw_pgam_Bu_new, mw_pgam_Bs_new)

# For VCM and int2D we only need the marginal B-splines
mw_Bu_new <- JOPS::bbase(mw_bins_u_new, 
                         xl = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                         xr = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                         nseg = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

mw_Bs_new <- JOPS::bbase(mw_bins_s_new, 
                         xl = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                         xr = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                         nseg = attr(r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

# Calculate hazards
#  Pgam
Eta_pgam_me <- me_pgam_B_new %*% r_mar_pgam_me$pgam_mar_me$alpha
Haz_pgam_me <- exp(Eta_pgam_me)
Haz_pgam_me <- matrix(Haz_pgam_me, length(me_bins_u_new), length(me_bins_s_new))

Eta_pgam_mw <- mw_pgam_B_new %*% r_mar_pgam_mw$pgam_mar_mw$alpha
Haz_pgam_mw <- exp(Eta_pgam_mw)
Haz_pgam_mw <- matrix(Haz_pgam_mw, length(mw_bins_u_new), length(mw_bins_s_new))

# VCM
eta_me <- me_Bs_new %*% r_mar_vcm_me$vcm_mar_me$Alpha[, 1]
vc_me  <- me_Bs_new %*% r_mar_vcm_me$vcm_mar_me$Alpha[, 2]
eta_mw <- mw_Bs_new %*% r_mar_vcm_mw$vcm_mar_mw$Alpha[, 1]
vc_mw  <- mw_Bs_new %*% r_mar_vcm_mw$vcm_mar_mw$Alpha[, 2]

HazVCM_me <- exp( matrix(rep(eta_me, length(me_bins_u_new)), 
                         byrow = T, nrow = length(me_bins_u_new), 
                         ncol = length(me_bins_s_new)) + 
                    outer(me_bins_u_new, as.vector(vc_me))  )
HazVCM_mw <- exp( matrix(rep(eta_mw, length(mw_bins_u_new)), 
                         byrow = T, nrow = length(mw_bins_u_new), 
                         ncol = length(mw_bins_s_new)) + 
                    outer(mw_bins_u_new, as.vector(vc_mw))  )

#  2D interaction
Haz2D_me <- exp(me_Bu_new %*% r_mar_nocov_18p_me$mod_no_cov$optimal_model$Alpha %*% t(me_Bs_new))
Haz2D_mw <- exp(mw_Bu_new %*% r_mar_nocov_18p_mw$mod_no_cov$optimal_model$Alpha %*% t(mw_Bs_new))


# ---- Figure comparing cross-sections of the hazard from the three models ----
legrange <- range(Haz_pgam_me, HazVCM_me, Haz2D_me,
                  Haz_pgam_mw, HazVCM_mw, Haz2D_mw)

mat <-  matrix(1:6, ncol=3, nrow = 2, byrow = T)
breaks <- seq(0, legrange[2]+.01, length=100)

par(mar= c(3,4,4,.5),
    oma = c(2, 0, 2,0),
    mgp = c(1.5,.5,0),
    cex.main = 1.5,
    cex.lab = 1.2,
    font.main = 1)
layout(mat)

grey7 <- grey.colors(n = 7)

matplot(x = mw_bins_s_new, y = t(Haz_pgam_mw), 
        type = "l", col = grey.colors(n = 7), #150), 
        lty = 1, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")
matplot(x = mw_bins_s_new, y = t(HazVCM_mw), 
        type = "l", col = grey.colors(n = 7), # = 150), 
        lty = 1, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")
matplot(x = mw_bins_s_new, y = t(Haz2D_mw), 
        type = "l", col = grey.colors(n = 7), #150), 
        lty = 1, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")


matplot(x = me_bins_s_new, y = t(Haz_pgam_me), 
        type = "l", 
        col = grey.colors(n = 7), # 150), 
        lty = 1, 
        ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")

# ---- legend goes here (bootom row, leftmost)
legend(x=0, y=0.33, ncol =2, bty="n", cex = 1.0,
       legend = paste(" u =", seq(20, 45, by=5)), lty=1, col=grey7)


matplot(x = me_bins_s_new, y = t(HazVCM_me), 
        type = "l", col = grey.colors(n =7), # = 150), 
        lty = 1, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")
matplot(x = me_bins_s_new, y = t(Haz2D_me), 
        type = "l", col = grey.colors(n = 7), #150), 
        lty = 1, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")


# Add common title for the first row
mtext("Men West Germany", 
      side = 3, line = -3.4, at = 0.52,
      outer = TRUE, cex = 1.25)
# Common title for second row
mtext("Men East Germany", 
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