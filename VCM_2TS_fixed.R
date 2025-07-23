# VCM_2TS fixed smoothing
VCM_2TS_fixed <- function(ExEv, binlist, nseg=10, bdeg=3, pord=2, kappa=1e-10,
                    itdetail = TRUE, lambdas){
  y      <- as.vector(ExEv$Y)
  expose <- as.vector(ExEv$R)
  u      <- rep(binlist$midu, ncol(ExEv$Y))
  s      <- rep(binlist$mids, each=nrow(ExEv$Y))
  
  m = length(y)
  sl <- min(binlist$bins_s)
  sr <- max(binlist$bins_s)
  
  # B-splines and penalty construction
  B = bbase(s, xl=sl, xr=sr, nseg = nseg, bdeg=bdeg)
  n = ncol(B)
  D = diff(diag(n), diff = pord)
  P = t(D) %*% D
  
  # Build VCM basis and block diag penalty matrix
  C = cbind(B, diag(u) %*% B)
  K = diag(2 * n) * kappa
  
  # Initialize
  eta = log((y + 0.1)/(expose + 0.2) )
  a = 0
 # lambdas = rep(1, 2)
  
  # Iterate for lambdas (Schall algorithm)
  for (it2 in 1:20) {
    Q = kronecker(diag(lambdas), P) + K
    
    # Fit VCM
    a1 = a
    for (it in 1:10) {
      mu = expose * exp(eta)
      z = y - mu + mu * eta
      W = diag(c(mu))
      S = t(C) %*% W %*% C
      anew = solve(S + Q, t(C) %*% z)
      da = max(abs(anew - a))
      if (da < 1e-4) break
      a = anew
      eta = C %*% a
    }
    
    # Compute effective dimensions and sums of squares of coefficients, aic, bic
    y[y == 0] <- 10 ^ -4
    mu_c <- mu
    mu_c[mu_c==0] <- 10 ^ -4
    dev <- 2 * sum(y * log(y / mu_c))
    G = solve(S + Q, S) 
    ssa = eds = rep(0, 2)
    g = diag(G)
    ed = sum(g)
    for (k in 1:2) {
      r = (k - 1) * n + (1:n)
      ssa[k] = sum(a[r] ^ 2)
      eds[k] = sum(g[r])
    }
    aic <- dev + 2 * ed
    n_obs <- sum(r > 0)
    bic <- dev + ed * log(n_obs)
    
    # lam = lambdas
    # lambdas = eds / ssa
    # dla = max(abs(log10(lam) - log10(lambdas)))
    # if (itdetail) cat(it2, dla, log10(lambdas), '\n' )
    
  }
  
  
  
  # Compute
  Alpha = matrix(a, n, 2)  # coefficients for log-baseline [, 1] and
  # for varying coefficient [, 2]
  vcm_results <- list(B=B, Alpha = Alpha, smoothpars = lambdas, eds = eds, 
                      ed = ed, aic = aic, bic = bic)
  return(vcm_results)
}
