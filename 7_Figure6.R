# ---- Figure 6 -----
# Hazard of separation over duration of the cohabitation
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
load("R_separation_vcm_me.Rda")
load("R_separation_vcm_mw.Rda")
load("R_int2d_separation_me.Rda")
load("R_int2d_separation_mw.Rda")
load("R_pgam_separation_me.Rda")
load("R_pgam_separation_mw.Rda")

# ----- Estimate hazard with denser grid -----
# First: make denser grid
# Second: re-evaluate B-splines
# Third: compute hazards

# Men East
bins_me <- r_pgam_sep_me$bins 
me_bins_u_new <- seq(20, 45, by=5)   
me_bins_s_new <- seq(min(bins_me$bins_s), max(bins_me$bins_s), by = .2)

# For Pgam model: these new intervals should be expanded to a grid
me_new_grid <- expand.grid(me_bins_u_new, me_bins_s_new)
names(me_new_grid) <- c("bins_u_new", "bins_s_new")

me_pgam_Bu_new <- JOPS::bbase(me_new_grid$bins_u_new, 
                              xl = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                              xr = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                              nseg = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

me_pgam_Bs_new <- JOPS::bbase(me_new_grid$bins_s_new, 
                              xl = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                              xr = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                              nseg = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

me_pgam_B_new <- cbind(me_pgam_Bu_new, me_pgam_Bs_new)

# For VCM and int2D we only need the marginal B-splines
me_Bu_new <- JOPS::bbase(me_bins_u_new, 
                         xl = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                         xr = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                         nseg = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

me_Bs_new <- JOPS::bbase(me_bins_s_new, 
                         xl = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                         xr = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                         nseg = attr(r_sep_nocov_18p_me$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

# Men West
bins_mw <- r_sep_pgam_mw$bins 
mw_bins_u_new <- seq(20, 45, by= 5)
mw_bins_s_new <- seq(min(bins_mw$bins_s), max(bins_mw$bins_s), by = .2)

# For Pgam model: these new intervals should be expanded to a grid
mw_new_grid <- expand.grid(mw_bins_u_new, mw_bins_s_new)
names(mw_new_grid) <- c("bins_u_new", "bins_s_new")

mw_pgam_Bu_new <- JOPS::bbase(mw_new_grid$bins_u_new, 
                              xl = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                              xr = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                              nseg = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

mw_pgam_Bs_new <- JOPS::bbase(mw_new_grid$bins_s_new, 
                              xl = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                              xr = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                              nseg = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

mw_pgam_B_new <- cbind(mw_pgam_Bu_new, mw_pgam_Bs_new)

# For VCM and int2D we only need the marginal B-splines
mw_Bu_new <- JOPS::bbase(mw_bins_u_new, 
                         xl = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "xl"),
                         xr = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "xr"),
                         nseg = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bu, "nseg"))

mw_Bs_new <- JOPS::bbase(mw_bins_s_new, 
                         xl = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "xl"),
                         xr = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "xr"),
                         nseg = attr(r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Bbases$Bs, "nseg"))

# Calculate hazards
#  Pgam
Eta_pgam_me <- me_pgam_B_new %*% r_pgam_sep_me$mod_pgam_sep_me$alpha
Eta_pgam_me_Matrix <- matrix(me_pgam_B_new %*% r_pgam_sep_me$mod_pgam_sep_me$alpha,
                             length(me_bins_u_new), length(me_bins_s_new))
Haz_pgam_me <- exp(Eta_pgam_me)
Haz_pgam_me <- matrix(Haz_pgam_me, length(me_bins_u_new), length(me_bins_s_new))

Eta_pgam_mw <- mw_pgam_B_new %*% r_sep_pgam_mw$pgam_sep_mw$alpha
Eta_pgam_mw_Matrix <- matrix(mw_pgam_B_new %*% r_sep_pgam_mw$pgam_sep_mw$alpha,
                             length(mw_bins_u_new), length(mw_bins_s_new))

Haz_pgam_mw <- exp(Eta_pgam_mw)
Haz_pgam_mw <- matrix(Haz_pgam_mw, length(mw_bins_u_new), length(mw_bins_s_new))

# VCM
eta_me <- me_Bs_new %*% r_sep_vcm_me$vcm_sep_me$Alpha[, 1]
vc_me  <- me_Bs_new %*% r_sep_vcm_me$vcm_sep_me$Alpha[, 2]

eta_mw <- mw_Bs_new %*% r_sep_vcm_mw$vcm_sep_mw$Alpha[, 1]
vc_mw  <- mw_Bs_new %*% r_sep_vcm_mw$vcm_sep_mw$Alpha[, 2]

HazVCM_me <- exp( matrix(rep(eta_me, length(me_bins_u_new)),
                         byrow = T, nrow = length(me_bins_u_new),
                         ncol = length(me_bins_s_new)) +
                    outer(me_bins_u_new, as.vector(vc_me))  )
HazVCM_mw <- exp( matrix(rep(eta_mw, length(mw_bins_u_new)),
                         byrow = T, nrow = length(mw_bins_u_new),
                         ncol = length(mw_bins_s_new)) +
                    outer(mw_bins_u_new, as.vector(vc_mw))  )

#  2D interaction
Haz2D_me <- exp(me_Bu_new %*% r_sep_nocov_18p_me$mod_no_cov$optimal_model$Alpha %*% t(me_Bs_new))
Haz2D_mw <- exp(mw_Bu_new %*% r_sep_nocov_18p_mw$mod_no_cov$optimal_model$Alpha %*% t(mw_Bs_new))


# ---- Figure ----
legrange <- range(Haz_pgam_me, Haz2D_me, 
                  Haz_pgam_mw, Haz2D_mw)

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
        type = "l", col = grey.colors(n = 7),  
        lty = 1, lwd = 1.3,
        ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")

matplot(x = mw_bins_s_new, y = t(Eta_pgam_mw_Matrix), 
        type = "l", col = grey.colors(n = 7),  
        lty = 1, lwd = 1.3,
        ylim = range(Eta_pgam_mw_Matrix),
        xlab = "Duration of the cohabitation",
        ylab = "Log-Hazard")


matplot(x = mw_bins_s_new, y = t(Haz2D_mw), 
        type = "l", col = grey.colors(n = 7),  
        lty = 1, lwd = 1.3,
        ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")

## -- East German Men --------------------
matplot(x = me_bins_s_new, y = (t(Haz_pgam_me)), 
        type = "l", 
        col = grey.colors(n = 7),  
        lty = 1,  lwd = 1.3,
        ylim = (c(0,legrange[2]+.01)),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")

# ---- legend goes here (bottom row, leftmost)
legend(x=0.0, y=0.30, ncol =2, bty="n", cex = 1.0,
       legend = paste(" u =", seq(20, 45, by=5)), lty=1, lwd = 1.3, col=grey7)

# --------  log-Hazard eta !!
matplot(x = me_bins_s_new, y = (t(Eta_pgam_me_Matrix)), 
        type = "l", 
        col = grey.colors(n = 7), 
        lty = 1, lwd = 1.3,
        ylim = range(Eta_pgam_mw_Matrix),
        xlab = "Duration of the cohabitation",
        ylab = "Log-Hazard")

matplot(x = me_bins_s_new, y = t(Haz2D_me), 
        type = "l", col = grey.colors(n = 7), 
        lty = 1, lwd = 1.3, ylim = c(0,legrange[2]+.01),
        xlab = "Duration of the cohabitation",
        ylab = "Hazard")

# Add common title for the first row
mtext("Men West Germany", 
      side = 3, line = -3.3, at = 0.52,
      outer = TRUE, cex = 1.15)
# Common title for second row
mtext("Men East Germany", 
      side = 3, line = -23.8, at = 0.52,
      outer = TRUE, cex = 1.15)

# Column titles
mtext("Hazard, log-additive", 
      side = 3, line = -0.7, at = 0.18, 
      outer = TRUE, cex = 1.2, font = 2)
mtext("Log-hazard, log-additive", 
      side = 3, line = -0.7, at = 0.53, 
      outer = TRUE, cex = 1.2, font = 2)
mtext("Hazard, 2D interaction", 
      side = 3, line = -0.7, at = 0.86, 
      outer = TRUE, cex = 1.2, font = 2)
