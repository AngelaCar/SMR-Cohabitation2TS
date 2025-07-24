# ---- Figure S5-----
# Hazard of marriage over age at entry and duration of the cohabitation
# comparison of 3 models for Women East and Women West

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
we_bins_u_new <- seq(min(bins_we$bins_u), max(bins_we$bins_u), by = .2)
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
ww_bins_u_new <- seq(min(bins_ww$bins_u), max(bins_ww$bins_u), by = .2)
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

# ---- Reproduce two panels of Fig4 ----
# Men East
# 
#windows(width=18, height=6)
#
legrange <- range(Haz_pgam_we, HazVCM_we, Haz2D_we,
                  Haz_pgam_ww, HazVCM_ww, Haz2D_ww)

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
image.plot(we_bins_u_new,
           we_bins_s_new, 
           Haz_pgam_we, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 25),
           main = "Women East - log-additive", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(we_bins_u_new,
        we_bins_s_new, 
        Haz_pgam_we,
        col = "white", labcex = 0.8,
        add=TRUE)
box()

# VCM
image.plot(we_bins_u_new,
           we_bins_s_new, 
           HazVCM_we, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 25),
           main = "Women East - varying coefficient", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(we_bins_u_new,
        we_bins_s_new, 
        HazVCM_we,
        col = "white", labcex = 0.8,
        add=TRUE)
box()

# 2D interaction
image.plot(we_bins_u_new,
           we_bins_s_new, 
           Haz2D_we, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 25),
           main = "Women East - 2D interaction", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(we_bins_u_new,
        we_bins_s_new, 
        Haz2D_we,
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
image.plot(ww_bins_u_new,
           ww_bins_s_new, 
           Haz_pgam_ww, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 25),
           main = "Women West - log-additive", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(ww_bins_u_new,
        ww_bins_s_new, 
        Haz_pgam_ww,
        col = "white", labcex = 0.8,
        add=TRUE)
box()

# VCM
image.plot(ww_bins_u_new,
           ww_bins_s_new, 
           HazVCM_ww, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 25),
           main = "Women West - varying coefficient", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(ww_bins_u_new,
        ww_bins_s_new, 
        HazVCM_ww,
        col = "white", labcex = 0.8,
        add=TRUE)
box()

# 2D interaction
image.plot(ww_bins_u_new,
           ww_bins_s_new, 
           Haz2D_ww, 
           col = pal_2ts(length(breaks)-1), 
           breaks = breaks,
           xlim = c(18, 47),
           ylim = c(0, 25),
           main = "Women West - 2D interaction", 
           xlab = "Age at entry into cohabitation",
           ylab = "Duration of the cohabitation")
contour(ww_bins_u_new,
        ww_bins_s_new, 
        Haz2D_ww,
        col = "white", labcex = 0.8,
        add=TRUE)
box()
