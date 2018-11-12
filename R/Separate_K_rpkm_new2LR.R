SeparateKRpkmNew2LR <- function(y1, y2, ROUNDs, Zcut, sig_cut, k0 = 2) {
  Cutted_id1 <- which(y1 < Zcut)
  Cutted_id2 <- which(y2 < Zcut)
  Cutted1 <- sum(y1 < Zcut)
  Cutted2 <- sum(y2 < Zcut)
  Cutted0 <- Cutted1 + Cutted2

  y1 <- y1[which(y1 >= Zcut)]
  y2 <- y2[which(y2 >= Zcut)]
  y0 <- c(y1, y2)

  ny1 <- Cutted1 + length(y1)
  ny2 <- Cutted2 + length(y2)

  sig1 <- sqrt(var(y1))
  sig2 <- sqrt(var(y2))

  k <- k0

  wl1_1 <- wl_1 <- rep(0, k)
  wl1_2 <- wl_2 <- rep(0, k)
  ul1_1 <- ul_1 <- rep(0, k)
  ul1_2 <- ul_2 <- rep(0, k)
  sigmal1_1 <- sigmal_1 <- rep(0, k)
  sigmal1_2 <- sigmal_2 <- rep(0, k)

  wl1_0 <- wl_0 <- rep(0, 2 * k - 1)
  ul1_0 <- ul_0 <- rep(0, 2 * k - 1)
  sigmal1_0 <- sigmal_0 <- rep(0, 2 * k - 1)

  n1 <- length(y1)
  ul0_1 <- c()
  for (i in 1:k) {
    tg_ic <- floor(i * n1/(k + 1))
    cc <- sort(y1)[tg_ic]
    if (tg_ic < 1) {
      cc <- sort(y1)[1] - 0.5
    }
    if (tg_ic > n1) {
      cc <- sort(y1)[n1] + 0.5
    }
    ul0_1 <- c(ul0_1, cc)
  }
  n2 <- length(y2)
  ul0_2 <- c()
  for (i in 1:k) {
    tg_ic <- floor(i * n2/(k + 1))
    cc <- sort(y2)[tg_ic]
    if (tg_ic < 1) {
      cc <- sort(y2)[1] - 0.5
    }
    if (tg_ic > n2) {
      cc <- sort(y2)[n2] + 0.5
    }
    ul0_2 <- c(ul0_2, cc)
  }

  wl0_1 <- rep(1/k, k)
  wl0_2 <- rep(1/k, k)
  sigmal0_1 <- rep(sig1, k)
  sigmal0_2 <- rep(sig2, k)

  wl0_0 <- rep(1/(2 * k - 1), 2 * k - 1)
  ul0_0 <- c(ul0_1, ul0_2[-1])
  sigmal0_0 <- rep(min(sig1, sig2), 2 * k - 1)

  wl_0 <- wl0_0
  wl_1 <- wl0_1
  wl_2 <- wl0_2
  ul_0 <- ul0_0
  ul_1 <- ul0_1
  ul_2 <- ul0_2
  sigmal_0 <- sigmal0_0
  sigmal_1 <- sigmal0_1
  sigmal_2 <- sigmal0_2

  table_theta_t1 <- matrix(0, length(y1), length(wl0_1))
  table_theta_t2 <- matrix(0, length(y2), length(wl0_2))
  table_theta_t0 <- matrix(0, c(length(y1) + length(y2)), length(wl0_0))

  w_table1 <- matrix(0, length(wl_1), ROUNDs)
  u_table1 <- matrix(0, length(ul_1), ROUNDs)
  sig_table1 <- matrix(0, length(sigmal_1), ROUNDs)

  w_table2 <- matrix(0, length(wl_2), ROUNDs)
  u_table2 <- matrix(0, length(ul_2), ROUNDs)
  sig_table2 <- matrix(0, length(sigmal_2), ROUNDs)

  w_table0 <- matrix(0, length(wl_0), ROUNDs)
  u_table0 <- matrix(0, length(ul_0), ROUNDs)
  sig_table0 <- matrix(0, length(sigmal_0), ROUNDs)


  ###########################
  rounds <- 1
  for (rounds in 1:ROUNDs) {

    w_table1[, rounds] <- wl_1
    u_table1[, rounds] <- ul_1
    sig_table1[, rounds] <- sigmal_1
    w_table2[, rounds] <- wl_2
    u_table2[, rounds] <- ul_2
    sig_table2[, rounds] <- sigmal_2
    w_table0[, rounds] <- wl_0
    u_table0[, rounds] <- ul_0
    sig_table0[, rounds] <- sigmal_0

    for (row in 1:nrow(table_theta_t1)) {
      all <- wl_1 * dnorm(y1[row], ul_1, sigmal_1)
      table_theta_t1[row, ] <- all/sum(all)
    }

    for (row in 1:nrow(table_theta_t2)) {
      all <- wl_2 * dnorm(y2[row], ul_2, sigmal_2)
      table_theta_t2[row, ] <- all/sum(all)
    }

    pZil1 <- Pi_Zj_Zcut_new(Zcut, ul_1, sigmal_1, wl_1)  ################################################################
    pZil2 <- Pi_Zj_Zcut_new(Zcut, ul_2, sigmal_2, wl_2)  ################################################################

    wl1 <- (apply(table_theta_t1, 2, sum) + pZil1 * Cutted1)/(nrow(table_theta_t1) + Cutted1)
    wl2 <- (apply(table_theta_t2, 2, sum) + pZil2 * Cutted2)/(nrow(table_theta_t2) + Cutted2)

    wl1 <- wl1/sum(wl1)
    wl2 <- wl2/sum(wl2)

    denoml <- (apply(table_theta_t1, 2, sum) + pZil1 * Cutted1)
    denom2 <- (apply(table_theta_t2, 2, sum) + pZil2 * Cutted2)

    for (id in 1:length(ul1_1)) {
      ul1_1[id] <- (sum(y1 * table_theta_t1[, id]) + (ul_1[id] - sigmal_1[id] * InverseMillsRatio(Zcut, ul_1, sigmal_1)[id]) * pZil1[id] * Cutted1)/denoml[id]
    }
    for (id in 1:length(ul1_2)) {
      ul1_2[id] <- (sum(y2 * table_theta_t2[, id]) + (ul_2[id] - sigmal_2[id] * InverseMillsRatio(Zcut, ul_2, sigmal_2)[id]) * pZil2[id] * Cutted2)/denom2[id]
    }

    for (id in 1:length(sigmal1_1)) {
      sigmal1_1[id] <- sqrt((sum((y1 - ul_1[id])^2 * table_theta_t1[, id]) + (sigmal_1[id])^2 * (1 - (Zcut - ul_1[id])/sigmal_1[id] * InverseMillsRatio(Zcut, ul_1, sigmal_1)[id]) * pZil1[id] * Cutted1)/denoml[id])
    }
    for (id in 1:length(sigmal1_2)) {
      sigmal1_2[id] <- sqrt((sum((y2 - ul_2[id])^2 * table_theta_t2[, id]) + (sigmal_2[id])^2 * (1 - (Zcut - ul_2[id])/sigmal_2[id] * InverseMillsRatio(Zcut, ul_2, sigmal_2)[id]) * pZil2[id] * Cutted2)/denom2[id])
    }

    #############################
    tg_R <- order(ul1_1)
    wl1 <- wl1[tg_R]
    ul1_1 <- ul1_1[tg_R]
    sigmal1_1 <- sigmal1_1[tg_R]

    tg_R <- order(ul1_2)
    wl2 <- wl2[tg_R]
    ul1_2 <- ul1_2[tg_R]
    sigmal1_2 <- sigmal1_2[tg_R]


    ccc1 <- wl1 * (ny1/(ny1 + ny2))
    ccc2 <- wl2 * (ny2/(ny1 + ny2))
    wad <- c(ccc1, ccc2[-1])
    wad[1] <- ccc1[1] + ccc2[1]
    wad <- wad/sum(wad)

    u_ad <- c(ul1_1, ul1_2[-1])
    u_ad[1] <- sum(c(ul1_1[1] * Cutted1/(Cutted1 + Cutted2), ul1_2[1] * Cutted2/(Cutted1 + Cutted2)))

    sig_ad <- c(sigmal1_1, sigmal1_2[-1])
    sig_ad[1] <- max(sigmal1_1[1], sigmal1_2[1])

    for (row in 1:nrow(table_theta_t0)) {
      all <- wad * dnorm(y0[row], u_ad, sig_ad)
      table_theta_t0[row, ] <- all/sum(all)
    }

    pZil0 <- Pi_Zj_Zcut_new(Zcut, u_ad, sig_ad, wad)  ################################################################

    wl0 <- (apply(table_theta_t0, 2, sum) + pZil0 * Cutted0)/(nrow(table_theta_t0) + Cutted0)
    wl0 <- wl0/sum(wl0)

    denom0 <- (apply(table_theta_t0, 2, sum) + pZil0 * Cutted0)

    for (id in 1:length(ul1_0)) {
      ul1_0[id] <- (sum(y0 * table_theta_t0[, id]) + (u_ad[id] - sig_ad[id] * InverseMillsRatio(Zcut, u_ad, sig_ad)[id]) * pZil0[id] * Cutted0)/denom0[id]
    }
    for (id in 1:length(sigmal1_0)) {
      sigmal1_0[id] <- sqrt((sum((y0 - u_ad[id])^2 * table_theta_t0[, id]) + (sig_ad[id])^2 * (1 - (Zcut - u_ad[id])/sig_ad[id] * InverseMillsRatio(Zcut, u_ad, sig_ad)[id]) * pZil0[id] * Cutted0)/denom0[id])
    }

    ul1_1[1] <- ul1_0[1]
    ul1_2[1] <- ul1_0[1]
    sigmal1_1[1] <- sigmal1_0[1]
    sigmal1_2[1] <- sigmal1_0[1]

    for (row in 1:nrow(table_theta_t1)) {
      all <- wl_1 * dnorm(y1[row], ul_1, sigmal_1)
      table_theta_t1[row, ] <- all/sum(all)
    }

    for (row in 1:nrow(table_theta_t2)) {
      all <- wl_2 * dnorm(y2[row], ul_2, sigmal_2)
      table_theta_t2[row, ] <- all/sum(all)
    }

    pZil1 <- Pi_Zj_Zcut_new(Zcut, ul_1, sigmal_1, wl_1)  ################################################################
    pZil2 <- Pi_Zj_Zcut_new(Zcut, ul_2, sigmal_2, wl_2)  ################################################################
    wl1 <- (apply(table_theta_t1, 2, sum) + pZil1 * Cutted1)/(nrow(table_theta_t1) + Cutted1)
    wl2 <- (apply(table_theta_t2, 2, sum) + pZil2 * Cutted2)/(nrow(table_theta_t2) + Cutted2)
    wl1 <- wl1/sum(wl1)
    wl2 <- wl2/sum(wl2)
    wl_0 <- wl0
    ul_0 <- u_ad
    sigmal_0 <- sig_ad

    #############################
    wl_1 <- wl1
    ul_1 <- ul1_1
    sigmal_1 <- sigmal1_1
    wl_2 <- wl2
    ul_2 <- ul1_2
    sigmal_2 <- sigmal1_2
    sigmal_1[1] <- min(sigmal_1[1], sig_cut)
    sigmal_2[1] <- min(sigmal_2[1], sig_cut)
    sigmal_1[2] <- min(sigmal_1[2], sig_cut)
    sigmal_2[2] <- min(sigmal_2[2], sig_cut)


    ul_1[1] <- min(ul_1[1], Zcut)
    ul_2[1] <- min(ul_2[1], Zcut)
    ul_1[2] <- max(ul_1[2], Zcut)
    ul_2[2] <- max(ul_2[2], Zcut)
  }
  rrr0 <- cbind(wl0, ul_0, sigmal_0)
  rrr1 <- cbind(wl_1, ul_1, sigmal_1)
  rrr2 <- cbind(wl_2, ul_2, sigmal_2)
  return(list(rrr0, rrr1, rrr2))
}


Separate_K_rpkm_new2LR <- function(y1, y2, ROUNDs, Zcut, sig_cut, k0 = 2) {
  Cutted_id1 <- which(y1 < Zcut)
  Cutted_id2 <- which(y2 < Zcut)
  Cutted1 <- sum(y1 < Zcut)
  Cutted2 <- sum(y2 < Zcut)
  Cutted0 <- Cutted1 + Cutted2

  y1 <- y1[which(y1 >= Zcut)]
  y2 <- y2[which(y2 >= Zcut)]
  y0 <- c(y1, y2)

  ny1 <- Cutted1 + length(y1)
  ny2 <- Cutted2 + length(y2)

  sig1 <- sqrt(var(y1))
  sig2 <- sqrt(var(y2))

  k <- k0

  wl1_1 <- wl_1 <- rep(0, k)
  wl1_2 <- wl_2 <- rep(0, k)
  ul1_1 <- ul_1 <- rep(0, k)
  ul1_2 <- ul_2 <- rep(0, k)
  sigmal1_1 <- sigmal_1 <- rep(0, k)
  sigmal1_2 <- sigmal_2 <- rep(0, k)

  wl1_0 <- wl_0 <- rep(0, 2 * k - 1)
  ul1_0 <- ul_0 <- rep(0, 2 * k - 1)
  sigmal1_0 <- sigmal_0 <- rep(0, 2 * k - 1)

  n1 <- length(y1)
  ul0_1 <- c()
  for (i in 1:k) {
    tg_ic <- floor(i * n1/(k + 1))
    cc <- sort(y1)[tg_ic]
    if (tg_ic < 1) {
      cc <- sort(y1)[1] - 0.5
    }
    if (tg_ic > n1) {
      cc <- sort(y1)[n1] + 0.5
    }
    ul0_1 <- c(ul0_1, cc)
  }
  n2 <- length(y2)
  ul0_2 <- c()
  for (i in 1:k) {
    tg_ic <- floor(i * n2/(k + 1))
    cc <- sort(y2)[tg_ic]
    if (tg_ic < 1) {
      cc <- sort(y2)[1] - 0.5
    }
    if (tg_ic > n2) {
      cc <- sort(y2)[n2] + 0.5
    }
    ul0_2 <- c(ul0_2, cc)
  }

  wl0_1 <- rep(1/k, k)
  wl0_2 <- rep(1/k, k)
  sigmal0_1 <- rep(sig1, k)
  sigmal0_2 <- rep(sig2, k)

  wl0_0 <- rep(1/(2 * k - 1), 2 * k - 1)
  ul0_0 <- c(ul0_1, ul0_2[-1])
  sigmal0_0 <- rep(min(sig1, sig2), 2 * k - 1)

  wl_0 <- wl0_0
  wl_1 <- wl0_1
  wl_2 <- wl0_2
  ul_0 <- ul0_0
  ul_1 <- ul0_1
  ul_2 <- ul0_2
  sigmal_0 <- sigmal0_0
  sigmal_1 <- sigmal0_1
  sigmal_2 <- sigmal0_2

  table_theta_t1 <- matrix(0, length(y1), length(wl0_1))
  table_theta_t2 <- matrix(0, length(y2), length(wl0_2))
  table_theta_t0 <- matrix(0, c(length(y1) + length(y2)), length(wl0_0))

  w_table1 <- matrix(0, length(wl_1), ROUNDs)
  u_table1 <- matrix(0, length(ul_1), ROUNDs)
  sig_table1 <- matrix(0, length(sigmal_1), ROUNDs)

  w_table2 <- matrix(0, length(wl_2), ROUNDs)
  u_table2 <- matrix(0, length(ul_2), ROUNDs)
  sig_table2 <- matrix(0, length(sigmal_2), ROUNDs)

  w_table0 <- matrix(0, length(wl_0), ROUNDs)
  u_table0 <- matrix(0, length(ul_0), ROUNDs)
  sig_table0 <- matrix(0, length(sigmal_0), ROUNDs)


  ###########################
  rounds <- 1
  for (rounds in 1:ROUNDs) {

    w_table1[, rounds] <- wl_1
    u_table1[, rounds] <- ul_1
    sig_table1[, rounds] <- sigmal_1
    w_table2[, rounds] <- wl_2
    u_table2[, rounds] <- ul_2
    sig_table2[, rounds] <- sigmal_2
    w_table0[, rounds] <- wl_0
    u_table0[, rounds] <- ul_0
    sig_table0[, rounds] <- sigmal_0

    for (i in 1:nrow(table_theta_t1)) {
      for (j in 1:ncol(table_theta_t1)) {
        table_theta_t1[i, j] <- sigma_f(wl_1, ul_1, sigmal_1, y1, i, j)
      }
    }
    for (i in 1:nrow(table_theta_t2)) {
      for (j in 1:ncol(table_theta_t2)) {
        table_theta_t2[i, j] <- sigma_f(wl_2, ul_2, sigmal_2, y2, i, j)
      }
    }

    pZil1 <- Pi_Zj_Zcut_new(Zcut, ul_1, sigmal_1, wl_1)  ################################################################
    pZil2 <- Pi_Zj_Zcut_new(Zcut, ul_2, sigmal_2, wl_2)  ################################################################

    wl1 <- (apply(table_theta_t1, 2, sum) + pZil1 * Cutted1)/(nrow(table_theta_t1) + Cutted1)
    wl2 <- (apply(table_theta_t2, 2, sum) + pZil2 * Cutted2)/(nrow(table_theta_t2) + Cutted2)

    wl1 <- wl1/sum(wl1)
    wl2 <- wl2/sum(wl2)

    denoml <- (apply(table_theta_t1, 2, sum) + pZil1 * Cutted1)
    denom2 <- (apply(table_theta_t2, 2, sum) + pZil2 * Cutted2)

    for (id in 1:length(ul1_1)) {
      ul1_1[id] <- (sum(y1 * table_theta_t1[, id]) + (ul_1[id] - sigmal_1[id] * InverseMillsRatio(Zcut, ul_1, sigmal_1)[id]) * pZil1[id] * Cutted1)/denoml[id]
    }
    for (id in 1:length(ul1_2)) {
      ul1_2[id] <- (sum(y2 * table_theta_t2[, id]) + (ul_2[id] - sigmal_2[id] * InverseMillsRatio(Zcut, ul_2, sigmal_2)[id]) * pZil2[id] * Cutted2)/denom2[id]
    }

    for (id in 1:length(sigmal1_1)) {
      sigmal1_1[id] <- sqrt((sum((y1 - ul_1[id])^2 * table_theta_t1[, id]) + (sigmal_1[id])^2 * (1 - (Zcut - ul_1[id])/sigmal_1[id] * InverseMillsRatio(Zcut, ul_1, sigmal_1)[id]) * pZil1[id] * Cutted1)/denoml[id])
    }
    for (id in 1:length(sigmal1_2)) {
      sigmal1_2[id] <- sqrt((sum((y2 - ul_2[id])^2 * table_theta_t2[, id]) + (sigmal_1[id])^2 * (1 - (Zcut - ul_2[id])/sigmal_2[id] * InverseMillsRatio(Zcut, ul_2, sigmal_2)[id]) * pZil2[id] * Cutted2)/denom2[id])
    }

    #############################
    tg_R <- order(ul1_1)
    wl1 <- wl1[tg_R]
    ul1_1 <- ul1_1[tg_R]
    sigmal1_1 <- sigmal1_1[tg_R]

    tg_R <- order(ul1_2)
    wl2 <- wl2[tg_R]
    ul1_2 <- ul1_2[tg_R]
    sigmal1_2 <- sigmal1_2[tg_R]


    ccc1 <- wl1 * (ny1/(ny1 + ny2))
    ccc2 <- wl2 * (ny2/(ny1 + ny2))
    wad <- c(ccc1, ccc2[-1])
    wad[1] <- ccc1[1] + ccc2[1]
    wad <- wad/sum(wad)

    u_ad <- c(ul1_1, ul1_2[-1])
    u_ad[1] <- sum(c(ul1_1[1] / n, ul1_2[1] / n))

    sig_ad <- c(sigmal1_1, sigmal1_2[-1])
    sig_ad[1] <- max(sigmal1_1[1], sigmal1_2[1])


    for (i in 1:nrow(table_theta_t0)) {
      for (j in 1:ncol(table_theta_t0)) {
        table_theta_t0[i, j] <- sigma_f(wad, u_ad, sig_ad, y0, i, j)
      }
    }

    pZil0 <- Pi_Zj_Zcut_new(Zcut, u_ad, sig_ad, wad)  ################################################################

    wl0 <- (apply(table_theta_t0, 2, sum) + pZil0 * Cutted0)/(nrow(table_theta_t0) + Cutted0)
    wl0 <- wl0/sum(wl0)

    denom0 <- (apply(table_theta_t0, 2, sum) + pZil0 * Cutted0)

    for (id in 1:length(ul1_0)) {
      ul1_0[id] <- (sum(y0 * table_theta_t0[, id]) + (u_ad[id] - sig_ad[id] * InverseMillsRatio(Zcut, u_ad, sig_ad)[id]) * pZil0[id] * Cutted0)/denom0[id]
    }
    for (id in 1:length(sigmal1_0)) {
      sigmal1_0[id] <- sqrt((sum((y0 - u_ad[id])^2 * table_theta_t0[, id]) + (sig_ad[id])^2 * (1 - (Zcut - u_ad[id])/sig_ad[id] * InverseMillsRatio(Zcut, u_ad, sig_ad)[id]) * pZil0[id] * Cutted0)/denom0[id])
    }

    ul1_1[1] <- ul1_0[1]
    ul1_2[1] <- ul1_0[1]
    sigmal1_1[1] <- sigmal1_0[1]
    sigmal1_2[1] <- sigmal1_0[1]

    for (i in 1:nrow(table_theta_t1)) {
      for (j in 1:ncol(table_theta_t1)) {
        table_theta_t1[i, j] <- sigma_f(wl_1, ul_1, sigmal_1, y1, i, j)
      }
    }
    for (i in 1:nrow(table_theta_t2)) {
      for (j in 1:ncol(table_theta_t2)) {
        table_theta_t2[i, j] <- sigma_f(wl_2, ul_2, sigmal_2, y2, i, j)
      }
    }

    pZil1 <- Pi_Zj_Zcut_new(Zcut, ul_1, sigmal_1, wl_1)  ################################################################
    pZil2 <- Pi_Zj_Zcut_new(Zcut, ul_2, sigmal_2, wl_2)  ################################################################
    wl1 <- (apply(table_theta_t1, 2, sum) + pZil1 * Cutted1)/(nrow(table_theta_t1) + Cutted1)
    wl2 <- (apply(table_theta_t2, 2, sum) + pZil2 * Cutted2)/(nrow(table_theta_t2) + Cutted2)
    wl1 <- wl1/sum(wl1)
    wl2 <- wl2/sum(wl2)
    wl_0 <- wl0
    ul_0 <- u_ad
    sigmal_0 <- sig_ad

    #############################
    wl_1 <- wl1
    ul_1 <- ul1_1
    sigmal_1 <- sigmal1_1
    wl_2 <- wl2
    ul_2 <- ul1_2
    sigmal_2 <- sigmal1_2
    # sigmal1_0[1]<-min(sigmal1_0[1],sig_cut)
    sigmal_1[1] <- min(sigmal_1[1], sig_cut)
    sigmal_2[1] <- min(sigmal_2[1], sig_cut)
    # sigmal1_0[2]<-min(sigmal1_0[2],sig_cut)
    sigmal_1[2] <- min(sigmal_1[2], sig_cut)
    sigmal_2[2] <- min(sigmal_2[2], sig_cut)


    # sigmal1_0[1]<-min(sigmal1_0[1],sig_cut)
    ul_1[1] <- min(ul_1[1], Zcut)
    ul_2[1] <- min(ul_2[1], Zcut)
    #ul_3[1] <- min(ul_3[1], Zcut)

    # sigmal1_0[2]<-min(sigmal1_0[2],sig_cut)
    ul_1[2] <- max(ul_1[2], Zcut)
    ul_2[2] <- max(ul_2[2], Zcut)
    #ul_3[2] <- max(ul_3[2], Zcut)

  }
  rrr0 <- cbind(wl0, ul_0, sigmal_0)
  rrr1 <- cbind(wl_1, ul_1, sigmal_1)
  rrr2 <- cbind(wl_2, ul_2, sigmal_2)
  return(list(rrr0, rrr1, rrr2))
}
