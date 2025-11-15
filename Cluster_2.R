## ============================================================
## C = 2 STUDY (10 parameter setups × 100 datasets each)
## - For each setup:
##   (A) Scatter montage of first 12 datasets with TRUE labels
##   (B) For each of 100 datasets:
##       * Fit AGM for C = 1:6, compute train/test joint RMSE
##       * Choose C* by visual-knee (Kneedle + L-method, fallback)
##       * (Optional) save elbow plot per dataset
##   (C) Histogram of the 100 chosen C*
## ============================================================

set.seed(2025)

## -------------------- KNOBS --------------------
n_per_dataset <- 300
n_datasets    <- 100         # per setup
Cs_grid       <- 1:6
train_frac    <- 0.7

# EM controls (kept light but stable)
em_max_iter   <- 120
em_restarts   <- 5
min_eig       <- 1e-4
min_weight    <- 2e-3

# Visual-knee controls
kneedle_span  <- 0.4         # 0.3–0.5 typical
fallback_eps  <- 0.01        # smallest C within (1+eps)*min(TEST)

# Plot/output controls
save_elbows   <- TRUE       # TRUE will write 100 elbow PNGs per setup (heavy)
scatter_panels_to_show <- 12 # how many datasets to show in scatter montage (first N)
panel_layout  <- c(3,4)      # rows × cols for scatter montage (match ^)
out_dir <- "AGM_C2_STUDY"
dir.create(out_dir, showWarnings = FALSE)

## -------------------- UTILITIES --------------------
rmvn2 <- function(n, mu, Sigma){
  Z <- matrix(rnorm(2*n), n, 2); L <- chol(Sigma); sweep(Z %*% L, 2, mu, `+`)
}
# mixture with labels
rmixture2_lbl <- function(n, pis, mus, Sigmas){
  C <- length(pis); z <- sample.int(C, n, TRUE, pis)
  Y <- matrix(NA_real_, n, 2)
  for(k in 1:C){ idx <- which(z==k); if(length(idx)) Y[idx,] <- rmvn2(length(idx), mus[k,], Sigmas[[k]]) }
  list(Y = Y, z = z)
}
log_dmvnorm2 <- function(y, mu, Sigma){
  L <- chol(Sigma); z <- t(backsolve(L, t(y)-mu, transpose=TRUE))
  -0.5*(2*log(2*pi) + 2*sum(log(diag(L))) + rowSums(z^2))
}
agm_em <- function(Y, C, max_iter=120, tol=1e-6, n_starts=5, seed=123){
  set.seed(seed)
  best <- NULL; best_ll <- -Inf
  for(s in 1:n_starts){
    km <- kmeans(Y, centers=C, nstart=1)
    n  <- nrow(Y)
    pis <- pmax(as.numeric(table(km$cluster))/n, min_weight); pis <- pis/sum(pis)
    mus <- km$centers
    Sigmas <- lapply(1:C, function(k){
      Yk <- Y[km$cluster==k,,drop=FALSE]
      S  <- if(nrow(Yk)<2) diag(apply(Y,2,stats::var)+1e-3) else stats::cov(Yk)
      ev <- eigen(S, symmetric=TRUE, only.values=TRUE)$values
      add <- max(0, min_eig - min(ev)); S + add*diag(2)
    })
    for(it in 1:max_iter){
      log_comp <- sapply(1:C, function(k) log(pis[k]) + log_dmvnorm2(Y, mus[k,], Sigmas[[k]]))
      m <- apply(log_comp, 1, max); log_den <- m + log(rowSums(exp(log_comp - m)))
      gamma <- exp(log_comp - log_den)
      Nk <- colSums(gamma) + 1e-12
      pis_new <- pmax(Nk/sum(Nk), min_weight); pis_new <- pis_new/sum(pis_new)
      mus_new <- t(gamma) %*% Y / Nk
      Sigmas_new <- lapply(1:C, function(k){
        Yc <- sweep(Y, 2, mus_new[k,], `-`); W <- sqrt(gamma[,k])
        S  <- crossprod(Yc*W)/Nk[k]
        ev <- eigen(S, symmetric=TRUE, only.values=TRUE)$values
        add <- max(0, min_eig - min(ev)); S + add*diag(2)
      })
      dpar <- sqrt(sum((pis_new-pis)^2) + sum((mus_new-mus)^2) +
                     sum(sapply(1:C, function(k) sum((Sigmas_new[[k]]-Sigmas[[k]])^2))))
      pis <- pis_new; mus <- mus_new; Sigmas <- Sigmas_new
      if(dpar < tol) break
    }
    log_comp <- sapply(1:C, function(k) log(pis[k]) + log_dmvnorm2(Y, mus[k,], Sigmas[[k]]))
    m <- apply(log_comp, 1, max); ll <- sum(m + log(rowSums(exp(log_comp - m))))
    if(ll > best_ll){ best <- list(pis=as.numeric(pis), mu=as.matrix(mus), Sigma=Sigmas, loglik=ll); best_ll <- ll }
  }
  best
}
resp <- function(Y, fit){
  C <- length(fit$pis)
  log_comp <- sapply(1:C, function(k) log(fit$pis[k]) + log_dmvnorm2(Y, fit$mu[k,], fit$Sigma[[k]]))
  m <- apply(log_comp, 1, max); den <- m + log(rowSums(exp(log_comp - m)))
  exp(log_comp - den)
}
predict_mean <- function(Y, fit) resp(Y, fit) %*% fit$mu
rmse_joint <- function(truth, pred) sqrt(mean(rowSums((truth - pred)^2)))

## Visual knee pickers
enforce_monotone_decreasing <- function(y) cummin(y)
kneedle_elbow <- function(Cs, vals, span = 0.4) {
  x <- as.numeric(Cs); y <- as.numeric(vals)
  xn <- (x - min(x)) / (max(x) - min(x))
  yn <- (y - min(y)) / max(1e-12, (max(y) - min(y)))
  yns <- as.numeric(predict(loess(yn ~ xn, span = span, degree = 1),
                            data.frame(xn = xn)))
  d <- xn - yns
  idx <- 2:(length(x)-1)
  x[ idx[ which.max(d[idx]) ] ]
}
L_method_elbow <- function(Cs, vals) {
  x <- as.numeric(Cs); y <- as.numeric(vals)
  K <- length(x); if (K < 3) return(x[which.min(y)])
  sse <- rep(Inf, K)
  for (k in 2:(K-1)) {
    ml <- lm(y[1:k] ~ x[1:k]); mr <- lm(y[k:K] ~ x[k:K])
    sse[k] <- sum(resid(ml)^2) + sum(resid(mr)^2)
  }
  x[which.min(sse)]
}
fallback_within <- function(Cs, vals, alpha = 0.01) {
  thr <- (1 + alpha) * min(vals)
  Cs[which(vals <= thr)][1]
}
choose_C_star <- function(Cs, test_vals, alpha = 0.01, span = 0.4) {
  y <- enforce_monotone_decreasing(test_vals)
  c1 <- kneedle_elbow(Cs, y, span = span)
  c2 <- L_method_elbow(Cs, y)
  cand <- sort(unique(na.omit(c(c1, c2))))
  cand <- cand[cand != min(Cs) & cand != max(Cs)]
  if (length(cand) > 0) return(cand[1])
  fallback_within(Cs, y, alpha = alpha)
}

## Rotate/shape covariance
rot_cov <- function(a,b,phi){
  R <- matrix(c(cos(phi),-sin(phi),sin(phi),cos(phi)),2,2)
  S <- diag(c(a^2,b^2))
  R %*% S %*% t(R)
}

## -------------------- TEN PARAMETRIC SETUPS (true C=2) --------------------
# We vary separation (delta), imbalance (pi), rotation (phi), and eccentricity.
S1 <- rot_cov(1.0,0.6, 0)
S2 <- rot_cov(1.0,0.6, 0)
setup_list <- list(
  list(name="S1_sep_medium_balanced",
       pis=c(0.5,0.5),
       mus=rbind(c( 3, 0), c(-3, 0)),
       Sigmas=list(S1,S2)),
  list(name="S2_sep_large_balanced",
       pis=c(0.5,0.5),
       mus=rbind(c( 4.5, 0), c(-4.5, 0)),
       Sigmas=list(S1,S2)),
  list(name="S3_sep_small_balanced",
       pis=c(0.5,0.5),
       mus=rbind(c( 2, 0), c(-2, 0)),
       Sigmas=list(S1,S2)),
  list(name="S4_imbalanced_70_30",
       pis=c(0.7,0.3),
       mus=rbind(c( 3.5, 0), c(-3.5, 0)),
       Sigmas=list(S1,S2)),
  list(name="S5_imbalanced_80_20",
       pis=c(0.8,0.2),
       mus=rbind(c( 4, 0), c(-4, 0)),
       Sigmas=list(S1,S2)),
  list(name="S6_rotated_cov_phi30",
       pis=c(0.5,0.5),
       mus=rbind(c( 3.2, 0.5), c(-3.2, -0.5)),
       Sigmas=list(rot_cov(1.1,0.5, pi/6), rot_cov(1.1,0.5,-pi/6))),
  list(name="S7_rotated_cov_phi45",
       pis=c(0.5,0.5),
       mus=rbind(c( 3.2, 0), c(-3.2, 0)),
       Sigmas=list(rot_cov(1.2,0.5, pi/4), rot_cov(1.2,0.5,-pi/4))),
  list(name="S8_high_variance",
       pis=c(0.5,0.5),
       mus=rbind(c( 3.8, 0), c(-3.8, 0)),
       Sigmas=list(rot_cov(1.6,1.0, 0), rot_cov(1.6,1.0, 0))),
  list(name="S9_eccentric_diff",
       pis=c(0.55,0.45),
       mus=rbind(c( 3.5, 0.5), c(-3.5, -0.5)),
       Sigmas=list(rot_cov(1.4,0.4, 0), rot_cov(0.8,0.6, 0))),
  list(name="S10_close_overlap",
       pis=c(0.5,0.5),
       mus=rbind(c( 2.6, 0.3), c(-2.6, -0.3)),
       Sigmas=list(rot_cov(1.2,0.8, 0), rot_cov(1.2,0.8, 0)))
)

## -------------------- MAIN LOOP OVER SETUPS --------------------
for (s in seq_along(setup_list)) {
  ps <- setup_list[[s]]
  sdir <- file.path(out_dir, ps$name)
  dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("\n=== %s ===\n", ps$name))
  
  ## (A) SCATTER MONTAGE (first `scatter_panels_to_show` datasets)
  scat_file <- file.path(sdir, "A_scatter_true_labels.png")
  png(scat_file, width=1600, height=1000)
  oldpar <- par(mfrow=panel_layout, mar=c(3.2,3.2,2.6,1.2), mgp=c(2.1,0.7,0))
  cols <- c("#1b9e77","#d95f02")
  to_show <- min(scatter_panels_to_show, n_datasets)
  for (d in 1:to_show) {
    set.seed(10000 + 100*s + d)
    sim <- rmixture2_lbl(n_per_dataset, ps$pis, ps$mus, ps$Sigmas)
    Y <- sim$Y; z <- sim$z
    plot(Y[,1], Y[,2], col = cols[z], pch = 16, cex = 0.7,
         xlab = "Y1", ylab = "Y2",
         main = sprintf("TRUE C=2 | ds=%d", d))
    grid()
    if (d == 1) legend("topright", legend=c("Comp 1","Comp 2"),
                       col = cols, pch=16, bty="n", cex=0.9)
  }
  par(oldpar); dev.off()
  
  ## (B) ELBOW + C* SELECTION FOR ALL 100 DATASETS
  Cstars <- integer(n_datasets)
  
  if (save_elbows) {
    elbows_dir <- file.path(sdir, "B_elbows")
    dir.create(elbows_dir, showWarnings = FALSE)
  }
  
  for (d in 1:n_datasets) {
    if (d %% 10 == 0) cat(sprintf(".. dataset %d/%d\n", d, n_datasets))
    set.seed(20000 + 100*s + d)
    sim <- rmixture2_lbl(n_per_dataset, ps$pis, ps$mus, ps$Sigmas)
    Y <- sim$Y
    
    idx_tr <- sample(seq_len(nrow(Y)), floor(train_frac*nrow(Y)))
    Ytr <- Y[idx_tr,,drop=FALSE]; Yte <- Y[-idx_tr,,drop=FALSE]
    
    rmse_tr <- rmse_te <- numeric(length(Cs_grid))
    for (i in seq_along(Cs_grid)) {
      C <- Cs_grid[i]
      fit <- agm_em(Ytr, C, max_iter=em_max_iter, n_starts=em_restarts,
                    seed = 777 + C + d*13 + s*1000)
      Yhat_tr <- predict_mean(Ytr, fit)
      Yhat_te <- predict_mean(Yte, fit)
      rmse_tr[i] <- rmse_joint(Ytr, Yhat_tr)
      rmse_te[i] <- rmse_joint(Yte, Yhat_te)
    }
    
    C_star <- choose_C_star(Cs_grid, rmse_te, alpha=fallback_eps, span=kneedle_span)
    Cstars[d] <- C_star
    
    if (save_elbows) {
      ylim <- range(c(rmse_tr, rmse_te))
      png(file.path(elbows_dir, sprintf("elbow_ds%03d.png", d)), width=900, height=650)
      plot(Cs_grid, rmse_tr, type="b", pch=19, ylim=ylim,
           xlab="C", ylab="Joint RMSE",
           main=sprintf("Elbow ds=%d | C*=%d", d, C_star))
      lines(Cs_grid, rmse_te, type="b", pch=17, lty=2)
      segments(Cs_grid, rmse_tr, Cs_grid, rmse_te, col=gray(0.85))
      abline(v = C_star, col="red", lty=3, lwd=2)
      legend("topright", c("Train RMSE","Test RMSE","Gap","Chosen C*"),
             pch=c(19,17,NA,NA), lty=c(1,2,1,3),
             col=c("black","black",gray(0.85),"red"), bty="n", cex=0.9)
      grid(); dev.off()
    }
  }
  
  # Save C* selections and histogram
  write.csv(data.frame(dataset=1:n_datasets, C_star=Cstars),
            file.path(sdir, "C_Cstar_selections.csv"), row.names=FALSE)
  
  png(file.path(sdir, "D_hist_Cstar.png"), width=900, height=600)
  freq <- table(factor(Cstars, levels = Cs_grid))
  bp <- barplot(freq, main=sprintf("%s: C* (visual knee; %d datasets)", ps$name, n_datasets),
                xlab="C*", ylab=sprintf("Count (out of %d)", n_datasets))
  text(bp, as.integer(freq), labels=as.integer(freq), pos=3, cex=0.9)
  abline(v = bp[which(Cs_grid == 2)], col="forestgreen", lty=2, lwd=2)
  mtext("True C = 2", side=3, line=0.2, col="forestgreen")
  dev.off()
  
  cat("  C* counts:\n")
  print(table(factor(Cstars, levels = Cs_grid)))
}

cat(sprintf("\nAll outputs saved in: %s\n", normalizePath(out_dir)))
