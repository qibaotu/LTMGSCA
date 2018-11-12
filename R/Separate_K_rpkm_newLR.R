SeparateKRpkmNewLR <- function(x, n, q, r, k = 2) {
  c <- colSums(x < q)
  c_sum <- sum(c)

  x_r <- apply(x, 2, function(x) x[which(x >= q)])

  x_all <- c(x_r, recursive = TRUE)

  mean0 <- matrix(nrow = k, ncol = ncol(x))
  sd0 <- matrix(nrow = k, ncol = ncol(x))

  mean <- matrix(nrow = k, ncol = ncol(x))
  for (col in 1:ncol(x)) {
    ni <- length(x_r[[col]])
    for (row in 1:k) {
      tg_ic <- floor(row * ni / (k + 1))
      cc <- sort(x_r[[col]])[tg_ic]
      if (tg_ic < 1) {
        cc <- sort(x_r[[col]])[1] - 0.5
      }
      if (tg_ic > ni) {
        cc <- sort(x_r[[col]])[ni] + 0.5
      }
      mean[row, col] <- cc
    }
  }

  p <- matrix(1 / k, nrow = k, ncol = ncol(x))
  sd <- matrix(sqrt(vapply(x_r, var, 0)), nrow = k, ncol = ncol(x), byrow = TRUE)
  t <- lapply(x_r, function(x) matrix(nrow = length(x), ncol = k))
  t0 <- matrix(nrow = sum(vapply(x_r, length, 0)), ncol = k + ncol(x) - 1)

  for (i in 1:n) {
    ccc <- matrix(nrow = k, ncol = ncol(x))
    wad <- rep(0, ncol(x) + 1)
    mean_all <- rep(0, ncol(x) + 1)
    sd_all <- rep(0, ncol(x) + 1)
    for (col in 1:ncol(x)) {
      t_u <- t(p[, col] * vapply(x_r[[col]], function(x) dnorm(x, mean[, col], sd[, col]), rep(0, k)))
      t[[col]] <- t_u / rowSums(t_u)
      pZil2 <- Pi_Zj_Zcut_new(q, mean[, col], sd[, col], p[, col])
      denom2 <- (colSums(t[[col]]) + pZil2 * c[[col]])
      im <- InverseMillsRatio(q, mean[, col], sd[, col])
      mean0[, col] <- colSums(crossprod(x_r[[col]], t[[col]]) + (mean[, col] - sd[, col] * im) * pZil2 * c[[col]]) / denom2
      sd0[, col] <- sqrt((colSums((x_r[[col]] - matrix(mean[, col], ncol = length(mean[, col]), nrow = length(x_r[[col]]), byrow = TRUE))^2 * t[[col]]) + (sd[, col])^2 * (1 - (q - mean[, col]) / sd[, col] * im) * pZil2 * c[[col]]) / denom2)
      wl2 <- denom2 / (nrow(t[[col]]) + c[[col]])

      # Reorder by mean0
      tg_R <- order(mean0[, col])
      ccc[, col] <- wl2[tg_R]
      mean0[, col] <- mean0[, col][tg_R]
      sd0[, col] <- sd0[, col][tg_R]

      wad[[col + 1]] <- ccc[-1, col]
      wad[1] <- wad[1] + ccc[1, col]
      mean_all[[col + 1]] <- mean0[, col][-1]
      mean_all[1] <- mean_all[1] + mean0[1, col] * c[[col]] / c_sum
      sd_all[[col + 1]] <- sd0[, col][-1]
      sd_all[1] <- max(sd_all[1], sd0[1, col])
    }

    for (row in 1:nrow(t0)) {
      t0_u <- wad * dnorm(x_all[row], mean_all, sd_all)
      t0[row, ] <- t0_u / sum(t0_u)
    }

    pZil0 <- Pi_Zj_Zcut_new(q, mean_all, sd_all, wad)
    denom0 <- colSums(t0) + pZil0 * c_sum
    p_all <- denom0 / (nrow(t0) + c_sum)

    im1 <- InverseMillsRatio(q, mean_all[1], sd_all[1])
    mean0[1, ] <- (sum(x_all * t0[, 1]) + (mean_all[1] - sd_all[1] * im1) * pZil0[1] * c_sum) / denom0[1]
    sd0[1, ] <- sqrt((sum((x_all - mean_all[1])^2 * t0[, 1]) + sd_all[1]^2 * (1 - (q - mean_all[1]) / sd_all[1] * im1) * pZil0[1] * c_sum) / denom0[1])

    for (col in 1:ncol(x)) {
      t_u <- t(p[, col] * vapply(x_r[[col]], function(x) dnorm(x, mean[, col], sd[, col]), c(0, 0)))
      t[[col]] <- t_u / rowSums(t_u)
      pZil <- Pi_Zj_Zcut_new(q, mean[, col], sd[, col], p[, col])
      p_u <- (colSums(t[[col]]) + pZil * c[[col]]) / (nrow(t[[col]]) + c[[col]])
      p[, col] <- p_u / sum(p_u)
    }

    mean <- mean0
    sd <- sd0

    for (col in 1:ncol(x)) {
      sd[1, col] <- min(sd[1, col], r)
      sd[1, col] <- min(sd[2, col], r)
      mean[1, col] <- min(mean[1, col], q)
      mean[2, col] <- max(mean[2, col], q)
    }
  }

  ret <- vector("list", ncol(x) + 1)

  ret[[0 + 1]] <- cbind(p_all, mean_all, sd_all)
  for (col in 1:ncol(x)) {
    ret[[col + 1]] <- cbind(p[, col], mean[, col], sd[, col])
  }
  return(ret)
}
