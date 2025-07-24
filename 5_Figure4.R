# ---- Figure 4 -----
# Hazard of marriage over age at entry and duration of the cohabitation
# comparison of 3 models for Men East and Men West

#
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# ------------------------------------------------
rm(list=ls())
# ------------------------------------------------


library(fields)
library(viridis)
library(TwoTimeScales)
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
me_bins_u_new <- seq(min(bins_me$bins_u), max(bins_me$bins_u), by = .2)
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
mw_bins_u_new <- seq(min(bins_mw$bins_u), max(bins_mw$bins_u), by = .2)
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

# ---- Reproduce two panels of Fig4 ----
# Men East
# 
#windows(width=18, height=6)
#
legrange <- range(Haz_pgam_me, HazVCM_me, Haz2D_me,
                  Haz_pgam_mw, HazVCM_mw, Haz2D_mw)

pal_2ts <- function(n, alpha){rev(plasma(n, alpha))}
mat <-  matrix(1:3, ncol=3, nrow = 1, byrow = T)
breaks <- seq(0, legrange[2]+.01, length=100)


par(mar = c(4, 4, 3.3, 1),
    oma = c(2.5, 0, 0, .5),
    mgp = c(2.5,.5,0),
    cex.main = 2.5,
    cex.lab = 2.2,
    font.main = 1)
layout(mat)

# pgam
image.plot(me_bins_u_new,
           me_bins_s_new, 
           Haz_pgam_me, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 26),
           main = "Men East - log-additive", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(me_bins_u_new,
        me_bins_s_new, 
        Haz_pgam_me,
        col = "white", labcex = 0.8,
        add=TRUE)
box()

# VCM
image.plot(me_bins_u_new,
           me_bins_s_new, 
           HazVCM_me, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 26),
           main = "Men East - varying coefficient", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(me_bins_u_new,
        me_bins_s_new, 
        HazVCM_me,
        col = "white", labcex = 0.8,
        add=TRUE)
box()

# 2D interaction
image.plot(me_bins_u_new,
           me_bins_s_new, 
           Haz2D_me, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 26),
           main = "Men East - 2D interaction", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(me_bins_u_new,
        me_bins_s_new, 
        Haz2D_me,
        col = "white", labcex = 0.8,
        add=TRUE)
box()

# Men West
# 
par(mar = c(4, 4, 3.3, 1),
    oma = c(2.5, 0, 0, .5),
    mgp = c(2.5,.5,0),
    cex.main = 2.5,
    cex.lab = 2.2,
    font.main = 1)
layout(mat)

# pgam
image.plot(mw_bins_u_new,
           mw_bins_s_new, 
           Haz_pgam_mw, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 25),
           main = "Men West - log-additive", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(mw_bins_u_new,
        mw_bins_s_new, 
        Haz_pgam_mw,
        col = "white", labcex = 0.8,
        add=TRUE)
box()

# VCM
image.plot(mw_bins_u_new,
           mw_bins_s_new, 
           HazVCM_mw, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 25),
           main = "Men West - varying coefficient", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(mw_bins_u_new,
        mw_bins_s_new, 
        HazVCM_mw,
        col = "white", labcex = 0.8,
        add=TRUE)
box()

# 2D interaction
image.plot(mw_bins_u_new,
           mw_bins_s_new, 
           Haz2D_mw, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 25),
           main = "Men West - 2D interaction", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(mw_bins_u_new,
        mw_bins_s_new, 
        Haz2D_mw,
        col = "white", labcex = 0.8,
        add=TRUE)
box()
