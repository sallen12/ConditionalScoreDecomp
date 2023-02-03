#' unconditional Brier score decomposition terms
#'
#' @param p vector of probability forecasts.
#' @param o vector of binary outcomes.
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
#'
#' bs_decomp(p = p, o = o, bins = 10)
#' bs_decomp(p = p_disc, o = o)
#' bs_decomp(p = p, o = o, bins = 10, method = "classical")
#' bs_decomp(p = p, o = o, method = "isotonic")
bs_decomp <- function(p, o, bins = NULL, method = "bias_corrected"){
  if(method == "isotonic"){
    terms <- bs_decomp_iso(p, o)
  }else if(method %in% c("classical", "bias_corrected")){
    if(is.null(bins)){
      bins <- length(unique(p))
      warning("Taking bins to be the number of unique values of p")
    }
    terms <- SpecsVerification::BrierDecomp(p, o, bins, bias.corrected = (method == "bias_corrected"))
    terms <- c(terms[1, c("UNC", "RES", "REL")], terms[1, "UNC"] - terms[1, "RES"] + terms[1, "REL"])
    names(terms) <- c("UNC", "RES", "REL", "TOT")
  }
  return(terms)
}

# function to obtain unconditional Brier score decomposition terms from isotonic regression
bs_decomp_iso <- function(p, o, bins, states = NULL){
  na_ind <- is.na(o)
  diag_obj <- reliabilitydiag::reliabilitydiag(x = p[!na_ind], y = o[!na_ind])
  terms <- c(summary(diag_obj)$uncertainty,
             summary(diag_obj)$discrimination,
             summary(diag_obj)$miscalibration,
             summary(diag_obj)$mean_score)
  names(terms) <- c("UNC", "RES", "REL", "TOT")
  return(terms)
}
