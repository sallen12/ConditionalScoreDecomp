#' conditional Brier score decomposition
#'
#' @param p vector of probability forecasts.
#' @param o vector of binary outcomes.
#' @param states vector of states on which to perform the conditional decomposition.
#' @param bins number of bins to use in the decomposition.
#' @param method string specifying which method should be used to estimate the decomposition terms.
#'
#' @return A named vector containing the decomposition terms.
#' @export
#'
#' @examples
#' p <- runif(100)
#' o <- rbinom(100, 1, 0.5)
#' p_disc <- sample(seq(0, 1, 0.1), 100, TRUE)
#' states = rep(c("G1", "G2"), each = 50)
#'
#' bs_decomp_cond(p = p, o = o, states = states, bins = 10)
#' bs_decomp_cond(p = p_disc, o = o, states = states)
#' bs_decomp_cond(p = p, o = o, states = states, bins = 10, method = "classical")
#' bs_decomp_cond(p = p, o = o, states = states, method = "isotonic")
bs_decomp_cond <- function(p, o, states, bins = NULL, method = "bias_corrected"){

  groups <- unique(states)
  n_j <- sapply(seq_along(groups), function(j) sum(states == groups[j]))
  N <- length(o)

  if(method == "isotonic"){

    terms_un <- bs_decomp_iso(p, o) # unconditional terms
    terms_cnd <- sapply(seq_along(groups), function(j)
      bs_decomp_iso(p[states == groups[j]], o[states == groups[j]])) # conditional terms

    unc_A <- sum((n_j/N)*terms_cnd['UNC', ])
    res_fA <- sum((n_j/N)*terms_cnd['RES', ])
    rel_fA <- sum((n_j/N)*terms_cnd['REL', ])

    res_A <- terms_un['UNC'] - unc_A
    res_Af <- rel_fA - terms_un['REL']

  }else if(method %in% c("classical", "bias_corrected")){

    if(is.null(bins)){
      bins <- length(unique(p))
      warning("Taking bins to be the number of unique values of p")
    }

    forecasts <- unique(p)
    obar <- mean(o)

    n_k <- sapply(seq_along(forecasts), function(k) sum(p == forecasts[k]))
    o_k <- sapply(seq_along(forecasts), function(k) sum(o[p == forecasts[k]]))
    obar_k <- sapply(seq_along(forecasts), function(k) mean(o[p == forecasts[k]]))

    o_j <- sapply(seq_along(groups), function(j) sum(o[states == groups[j]]))
    obar_j <- sapply(seq_along(groups), function(j) mean(o[states == groups[j]]))

    n_kj <- sapply(seq_along(forecasts), function(k)
      sapply(seq_along(groups), function(j) sum(p == forecasts[k] & states == groups[j])))
    o_kj <- sapply(seq_along(forecasts), function(k)
      sapply(seq_along(groups), function(j) sum(o[p == forecasts[k] & states == groups[j]])))
    obar_kj <- sapply(seq_along(forecasts), function(k)
      sapply(seq_along(groups), function(j) mean(o[p == forecasts[k] & states == groups[j]])))

    unc_A <- sum((n_j/N)*obar_j*(1 - obar_j))
    res_A <- sum((n_j/N)*((obar_j - obar)^2))
    res_fA <- sum((n_kj/N)*((obar_j - obar_kj)^2), na.rm = T)
    res_Af <- sum(t(n_kj/N)*((obar_k - t(obar_kj))^2), na.rm = T)
    rel_fA <- sum(t(n_kj/N)*((forecasts - t(obar_kj))^2), na.rm = T)
    tot <- unc_A - res_fA + rel_fA

    if(method == "bias_corrected"){
      bias1 <- sum(n_j*obar_j*(1 - obar_j)/(n_j - 1))/N
      bias2 <- obar*(1 - obar)/(N - 1)
      bias3 <- sum(n_kj*obar_kj*(1 - obar_kj)/(n_kj - 1), na.rm = T)/N
      bias4 <- sum(n_k*obar_k*(1 - obar_k)/(n_k - 1))/N

      unc_A <- unc_A + bias1
      res_A <- res_A - bias1 + bias2
      res_fA <- res_fA + bias1 - bias3
      res_Af <- res_Af - bias3 + bias4
      rel_fA <- rel_fA - bias3
    }

  }

  tot <- unc_A - res_fA + rel_fA
  terms <- c(unc_A, res_A, res_fA, res_Af, rel_fA, tot)
  names(terms) <- c("UNC_A", "RES_A", "RES_F|A", "RES_A|F", "REL_F|A", "TOT")

  return(terms)

}
