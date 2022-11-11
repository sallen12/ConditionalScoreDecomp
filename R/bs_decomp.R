# function to obtain unconditional Brier score decomposition terms
bs_decomp <- function(p, o, method = "bias_corrected"){

  p <- p[!is.na(o)]
  o <- o[!is.na(o)]

  if(method == "isotonic"){

    terms <- bs_decomp_iso(p, o)

  }else if(method %in% c("classical", "bias_corrected")){

    forecasts <- unique(p)
    N <- length(o)
    obar <- mean(o)

    n_k <- sapply(seq_along(forecasts), function(k) sum(p == forecasts[k]))
    o_k <- sapply(seq_along(forecasts), function(k) sum(o[p == forecasts[k]]))
    obar_k <- sapply(seq_along(forecasts), function(k) mean(o[p == forecasts[k]]))

    unc <- obar*(1 - obar)
    res <- sum((n_k/N)*((obar_k - obar)^2))
    rel <- sum((n_k/N)*((obar_k - forecasts)^2))

    if(method == "bias_corrected"){
      unc <- unc + unc/(N - 1)
      res <- res + unc/(N - 1) - sum(n_k*obar_k*(1 - obar_k)/(n_k - 1))/N
      rel <- rel - sum(n_k*obar_k*(1 - obar_k)/(n_k - 1))/N
    }

    tot <- unc - res + rel
    terms <- c(unc, res, rel, tot)
    names(terms) <- c("UNC", "RES", "REL", "TOT")

  }

  return(terms)

}

# function to obtain unconditional Brier score decomposition terms from isotonic regression
bs_decomp_iso <- function(p, o, states = NULL){
  na_ind <- is.na(o)

  diag_obj <- reliabilitydiag::reliabilitydiag(x = p[!na_ind], y = o[!na_ind])
  terms <- c(summary(diag_obj)$uncertainty,
             summary(diag_obj)$discrimination,
             summary(diag_obj)$miscalibration,
             summary(diag_obj)$mean_score)
  names(terms) <- c("UNC", "RES", "REL", "TOT")
  return(terms)
}
