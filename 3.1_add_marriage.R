# ---- Transitions from cohabitation to marriage ----
# Log-additive models

#---- Load the necessary libraries and functions ----

library(data.table)
library(TwoTimeScales)
source("iwls_1d.R")

# ---- Load the data ----
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

# ---- P-GAM ----
#* Men east ----
du <- ds <- .5 # half-year
data2d_me <- prepare_data(
  u = MenEast$agecoh_y,
  s_in = MenEast$dur_entry,
  s_out = MenEast$durcoh_y,
  events = as.numeric(MenEast$event_type == "marriage"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47
)

r <- c(data2d_me$bindata$R)
y <- c(data2d_me$bindata$Y)
midu <- data2d_me$bins$midu
mids <- data2d_me$bins$mids
grid <- expand.grid(midu, mids)
data.grid <- cbind(r, y, grid)
names(data.grid) <- c("r", "y", "midu", "mids")

nu <- length(midu)
ns <- length(mids)
## Parameters of the B-splines and their calculation

# number of segments on the two axes -
# it depends on the how many mid-points
nsu <- 11
nss <- 10
deg <- 3 # degree of the B-splines
pord <- 2 # order of the penalty

# Construct the marginal bases
Bu <- JOPS::bbase(data.grid$midu, xl = 18, xr = 47, nseg = nsu, bdeg = deg)
nbu <- dim(Bu)[2] # number of B-splines over u
Bs <- JOPS::bbase(data.grid$mids, xl = 0, xr = 26, nseg = nss, bdeg = deg)
nbs <- dim(Bs)[2] # number of B-splines over s
B <- cbind(Bu, Bs)
nb <- ncol(B)

Du <- diff(diag(nbu), diff = 2)
DutDu <- t(Du) %*% Du
Ds <- diff(diag(nbs), diff = 2)
DstDs <- t(Ds) %*% Ds

lru <- seq(-1, 3, by = .5)
lrs <- lru
Aic <- matrix(NA, length(lru), length(lrs))

for (i in 1:length(lru)) {
  for (j in 1:length(lrs)) {
    ru <- 10^lru[i]
    rs <- 10^lrs[j]
    Pu <- ru * DutDu
    Ps <- rs * DstDs
    P <- matrix(0, nb, nb)
    P[1:nbu, 1:nbu] <- Pu
    P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
    PpR <- P + 1e-6 * diag(nb)

    mod <- iwls_1d(r, y, Bs = B, P = PpR)
    Aic[i, j] <- mod$aic
  }
}

# optimal model
wlr <- which(Aic == min(Aic), arr.ind = T)
optim_lr <- c(lru[wlr[1]], lrs[wlr[2]])
optim_r <- 10^optim_lr
Pu <- optim_r[1] * DutDu
Ps <- optim_r[2] * DstDs
P <- matrix(0, nb, nb)
P[1:nbu, 1:nbu] <- Pu
P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
PpR <- P + 1e-6 * diag(nb)

mod_opt <- iwls_1d(r, y, Bs = B, P = PpR)
pgam_mar_me <- mod_opt

mu <- r * exp(B %*% mod_opt$alpha)
w <- c(mu)

# Calculate AIC
BWB <- t(B) %*% (w * B)
BWBpP <- BWB + PpR
K <- solve(BWBpP) %*% BWB

# compute ED
K11 <- K[1:nbu, 1:nbu]
K22 <- K[(nbu + 1):(nbu + nbs), (nbu + 1):(nbu + nbs)]
EDu <- sum(diag(K11))
EDs <- sum(diag(K22))

r_mar_pgam_me <- list(
  "pgam_mar_me" = mod_opt,
  ED_rho = list(
    lrhou = optim_lr[1],
    lrhos = optim_lr[2],
    EDu = EDu,
    EDs = EDs
  ),
  AIC = mod_opt$aic,
  Bsplines = list(
    Bu = Bu,
    Bs = Bs,
    B = B
  ),
  bins = data2d_me$bins
)
save(r_mar_pgam_me, file = "R_pgam_marriage_me.Rda")

#* Men West ----
du <- ds <- .5 # half-year
data2d_mw <- prepare_data(
  u = MenWest$agecoh_y,
  s_in = MenWest$dur_entry,
  s_out = MenWest$durcoh_y,
  events = as.numeric(MenWest$event_type == "marriage"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47,
  min_s = 0, max_s = 25
)

# Additive model
r <- c(data2d_mw$bindata$R)
y <- c(data2d_mw$bindata$Y)
midu <- data2d_mw$bins$midu
mids <- data2d_mw$bins$mids
grid <- expand.grid(midu, mids)
data.grid <- cbind(r, y, grid)
names(data.grid) <- c("r", "y", "midu", "mids")

nu <- length(midu)
ns <- length(mids)
## Parameters of the B-splines and their calculation

# number of segments on the two axes -
# it depends on the how many mid-points
nsu <- 11
nss <- 9
deg <- 3 # degree of the B-splines
pord <- 2 # order of the penalty

# Construct the marginal bases
Bu <- JOPS::bbase(data.grid$midu, xl = 18, xr = 47, nseg = nsu, bdeg = deg)
nbu <- dim(Bu)[2] # number of B-splines over u
Bs <- JOPS::bbase(data.grid$mids, xl = 0, xr = 25, nseg = nss, bdeg = deg)
nbs <- dim(Bs)[2] # number of B-splines over s
B <- cbind(Bu, Bs)
nb <- ncol(B)

Du <- diff(diag(nbu), diff = 2)
DutDu <- t(Du) %*% Du
Ds <- diff(diag(nbs), diff = 2)
DstDs <- t(Ds) %*% Ds

lru <- seq(-1, 3, by = .5)
lrs <- lru
Aic <- matrix(NA, length(lru), length(lrs))

for (i in 1:length(lru)) {
  for (j in 1:length(lrs)) {
    ru <- 10^lru[i]
    rs <- 10^lrs[j]
    Pu <- ru * DutDu
    Ps <- rs * DstDs
    P <- matrix(0, nb, nb)
    P[1:nbu, 1:nbu] <- Pu
    P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
    PpR <- P + 1e-6 * diag(nb)

    mod <- iwls_1d(r, y, Bs = B, P = PpR)
    Aic[i, j] <- mod$aic
  }
}

# optimal model
wlr <- which(Aic == min(Aic), arr.ind = T)
optim_lr <- c(lru[wlr[1]], lrs[wlr[2]])
optim_r <- 10^optim_lr
Pu <- optim_r[1] * DutDu
Ps <- optim_r[2] * DstDs
P <- matrix(0, nb, nb)
P[1:nbu, 1:nbu] <- Pu
P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
PpR <- P + 1e-6 * diag(nb)

mod_opt <- iwls_1d(r, y, Bs = B, P = PpR)
pgam_mar_mw <- mod_opt

mu <- r * exp(B %*% mod_opt$alpha)
w <- c(mu)

# Calculate AIC
BWB <- t(B) %*% (w * B)
BWBpP <- BWB + PpR
K <- solve(BWBpP) %*% BWB

# compute ED
K11 <- K[1:nbu, 1:nbu]
K22 <- K[(nbu + 1):(nbu + nbs), (nbu + 1):(nbu + nbs)]
EDu <- sum(diag(K11))
EDs <- sum(diag(K22))

r_mar_pgam_mw <- list(
  "pgam_mar_mw" = mod_opt,
  ED_rho = list(
    lrhou = optim_lr[1],
    lrhos = optim_lr[2],
    EDu = EDu,
    EDs = EDs
  ),
  AIC = mod_opt$aic,
  Bsplines = list(
    Bu = Bu,
    Bs = Bs,
    B = B
  ),
  bins = data2d_mw$bins
)
save(r_mar_pgam_mw, file = "R_pgam_marriage_mw.Rda")

#* Women East ----
du <- ds <- .5 # half-year
data2d_we <- prepare_data(
  u = WomEast$agecoh_y,
  s_in = WomEast$dur_entry,
  s_out = WomEast$durcoh_y,
  events = as.numeric(WomEast$event_type == "marriage"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47
)

# Additive model
r <- c(data2d_we$bindata$R)
y <- c(data2d_we$bindata$Y)
midu <- data2d_we$bins$midu
mids <- data2d_we$bins$mids
grid <- expand.grid(midu, mids)
data.grid <- cbind(r, y, grid)
names(data.grid) <- c("r", "y", "midu", "mids")

nu <- length(midu)
ns <- length(mids)
## Parameters of the B-splines and their calculation

# number of segments on the two axes -
# it depends on the how many mid-points
nsu <- 11
nss <- 10
deg <- 3 # degree of the B-splines
pord <- 2 # order of the penalty

# Construct the marginal bases
Bu <- JOPS::bbase(data.grid$midu, xl = 18, xr = 47, nseg = nsu, bdeg = deg)
nbu <- dim(Bu)[2] # number of B-splines over u
Bs <- JOPS::bbase(data.grid$mids, xl = 0, xr = 25, nseg = nss, bdeg = deg)
nbs <- dim(Bs)[2] # number of B-splines over s
B <- cbind(Bu, Bs)
nb <- ncol(B)

Du <- diff(diag(nbu), diff = 2)
DutDu <- t(Du) %*% Du
Ds <- diff(diag(nbs), diff = 2)
DstDs <- t(Ds) %*% Ds

lru <- seq(-1, 3, by = .5)
lrs <- lru
Aic <- matrix(NA, length(lru), length(lrs))

for (i in 1:length(lru)) {
  for (j in 1:length(lrs)) {
    ru <- 10^lru[i]
    rs <- 10^lrs[j]
    Pu <- ru * DutDu
    Ps <- rs * DstDs
    P <- matrix(0, nb, nb)
    P[1:nbu, 1:nbu] <- Pu
    P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
    PpR <- P + 1e-6 * diag(nb)

    mod <- iwls_1d(r, y, Bs = B, P = PpR)
    Aic[i, j] <- mod$aic
  }
}

# optimal model
wlr <- which(Aic == min(Aic), arr.ind = T)
optim_lr <- c(lru[wlr[1]], lrs[wlr[2]])
optim_r <- 10^optim_lr
Pu <- optim_r[1] * DutDu
Ps <- optim_r[2] * DstDs
P <- matrix(0, nb, nb)
P[1:nbu, 1:nbu] <- Pu
P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
PpR <- P + 1e-6 * diag(nb)

mod_opt <- iwls_1d(r, y, Bs = B, P = PpR)
pgam_sep_ww <- mod_opt

mu <- r * exp(B %*% mod_opt$alpha)
w <- c(mu)

# Calculate AIC
BWB <- t(B) %*% (w * B)
BWBpP <- BWB + PpR
K <- solve(BWBpP) %*% BWB

# compute ED
K11 <- K[1:nbu, 1:nbu]
K22 <- K[(nbu + 1):(nbu + nbs), (nbu + 1):(nbu + nbs)]
EDu <- sum(diag(K11))
EDs <- sum(diag(K22))

r_pgam_marr_we <- list(
  "mod_pgam_marr_we" = mod_opt,
  ED_rho = list(
    lrhou = optim_lr[1],
    lrhos = optim_lr[2],
    EDu = EDu,
    EDs = EDs
  ),
  AIC = mod_opt$aic,
  Bsplines = list(
    Bu = Bu,
    Bs = Bs,
    B = B
  ),
  bins = data2d_we$bins
)
save(r_pgam_marr_we, file = "R_pgam_marriage_we.Rda")

#* Women West ----
du <- ds <- .5 # half-year
data2d_ww <- prepare_data(
  u = WomWest$agecoh_y,
  s_in = WomWest$dur_entry,
  s_out = WomWest$durcoh_y,
  events = as.numeric(WomWest$event_type == "marriage"),
  ds = ds,
  du = du,
  min_u = 18, max_u = 47,
  min_s = 0, max_s = 27
)

# Additive model
r <- c(data2d_ww$bindata$R)
y <- c(data2d_ww$bindata$Y)
midu <- data2d_ww$bins$midu
mids <- data2d_ww$bins$mids
grid <- expand.grid(midu, mids)
data.grid <- cbind(r, y, grid)
names(data.grid) <- c("r", "y", "midu", "mids")

nu <- length(midu)
ns <- length(mids)
## Parameters of the B-splines and their calculation

# number of segments on the two axes -
# it depends on the how many mid-points
nsu <- 11
nss <- 10
deg <- 3 # degree of the B-splines
pord <- 2 # order of the penalty

# Construct the marginal bases
Bu <- JOPS::bbase(data.grid$midu, xl = 18, xr = 47, nseg = nsu, bdeg = deg)
nbu <- dim(Bu)[2] # number of B-splines over u
Bs <- JOPS::bbase(data.grid$mids, xl = 0, xr = 25, nseg = nss, bdeg = deg)
nbs <- dim(Bs)[2] # number of B-splines over s
B <- cbind(Bu, Bs)
nb <- ncol(B)

Du <- diff(diag(nbu), diff = 2)
DutDu <- t(Du) %*% Du
Ds <- diff(diag(nbs), diff = 2)
DstDs <- t(Ds) %*% Ds

lru <- seq(-1, 3, by = .5)
lrs <- lru
Aic <- matrix(NA, length(lru), length(lrs))

for (i in 1:length(lru)) {
  for (j in 1:length(lrs)) {
    ru <- 10^lru[i]
    rs <- 10^lrs[j]
    Pu <- ru * DutDu
    Ps <- rs * DstDs
    P <- matrix(0, nb, nb)
    P[1:nbu, 1:nbu] <- Pu
    P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
    PpR <- P + 1e-6 * diag(nb)

    mod <- iwls_1d(r, y, Bs = B, P = PpR)
    Aic[i, j] <- mod$aic
  }
}

# optimal model
wlr <- which(Aic == min(Aic), arr.ind = T)
optim_lr <- c(lru[wlr[1]], lrs[wlr[2]])
optim_r <- 10^optim_lr
Pu <- optim_r[1] * DutDu
Ps <- optim_r[2] * DstDs
P <- matrix(0, nb, nb)
P[1:nbu, 1:nbu] <- Pu
P[(nbu + 1):(nbs + nbu), (nbu + 1):(nbs + nbu)] <- Ps
PpR <- P + 1e-6 * diag(nb)

mod_opt <- iwls_1d(r, y, Bs = B, P = PpR)
pgam_sep_ww <- mod_opt

mu <- r * exp(B %*% mod_opt$alpha)
w <- c(mu)

# Calculate AIC
BWB <- t(B) %*% (w * B)
BWBpP <- BWB + PpR
K <- solve(BWBpP) %*% BWB

# compute ED
K11 <- K[1:nbu, 1:nbu]
K22 <- K[(nbu + 1):(nbu + nbs), (nbu + 1):(nbu + nbs)]
EDu <- sum(diag(K11))
EDs <- sum(diag(K22))

r_pgam_marr_ww <- list(
  "mod_pgam_marr_ww" = mod_opt,
  ED_rho = list(
    lrhou = optim_lr[1],
    lrhos = optim_lr[2],
    EDu = EDu,
    EDs = EDs
  ),
  AIC = mod_opt$aic,
  Bsplines = list(
    Bu = Bu,
    Bs = Bs,
    B = B
  ),
  bins = data2d_ww$bins
)
save(r_pgam_marr_ww, file = "R_pgam_marriage_ww.Rda")
