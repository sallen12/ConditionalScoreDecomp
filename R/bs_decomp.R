#' Decompositions of the Brier score
#'
#' Calculate decompositions of the Brier score into uncertainty, resolution, and
#' reliability components. A conditional decomposition can also be performed that
#' calculates these terms conditional on some events, or states having occurred.
#'
#' @param o vector of binary outcomes.
#' @param p vector of probability forecasts.
#' @param states vector of states on which to perform the conditional decomposition.
#' @param bins integer; number of bins to use in the decomposition.
#' @param method string; the method to be used to estimate the decomposition terms.
#'
#' @details
#' The vectors \code{p}, \code{o}, and \code{states} (if used) should have the
#' same length.
#'
#' \code{o} is a numeric vector of values that are either 0 or 1, representing
#' the observed values.
#' \code{p} is a numeric vector of values between 0 and 1, representing the
#' corresponding probability forecasts.
#' \code{states} is a character vector, with \code{states[i]} corresponding to
#' the state that occurs when \code{p[i]} is forecast and \code{o[i]} is observed.
#' The number of unique states should be small compared to the number of observations.
#'
#' The \code{bins} argument specifies how many bins will be used when calculating
#' the decomposition. The default is the number of unique elements in \code{p}.
#' A warning is returned if this is large compared to the number of observations.
#' This argument is not required when the decomposition is calculated using
#' isotonic regression (\code{method = "isotonic"}).
#'
#' The \code{method} argument takes three possible options: \code{method = "classical"},
#' performs the classical decomposition of the Brier score that was proposed by
#' Murphy (1973); \code{method = "bias-corrected"} performs the bias-corrected
#' decomposition proposed by Ferro and Fricker (2012), which is more appropriate
#' for small sample sizes; \code{method = "isotonic"} (the default) performs the
#' decomposition based on isotonic regression proposed by Dimitriadis et al. (2021).
#' See references below. Note that the bias-corrected approach may return negative
#' estimates of the three terms.
#'
#' Calculation of the Brier score decompositions using isotonic regression
#' leverages the \pkg{reliabilitydiag} package.
#'
#' @return A named vector containing the decomposition terms.
#'
#' @author Sam Allen
#'
#' @references
#' \emph{Decomposition of the Brier score}
#'
#' Murphy, A.H. (1973): `A new vector partition of the probability score'. \emph{Journal of Applied Meteorology and Climatology} 12, 595-600. \doi{10.1175/1520-0450(1973)012<0595:ANVPOT>2.0.CO;2}
#'
#' \emph{Bias-corrected decomposition of the Brier score}
#'
#' Ferro, C.A.T. and T.E. Fricker (2012): `A bias‐corrected decomposition of the Brier score. \emph{Quarterly Journal of the Royal Meteorological Society} 138, 1954-1960. \doi{10.1002/qj.1924}
#'
#' \emph{Isotonic regression-based decomposition of the Brier score}
#'
#' Dimitriadis, T., Gneiting, T. and A.I. Jordan (2021): `Stable reliability diagrams for probabilistic classifiers'. \emph{Proceedings of the National Academy of Sciences} 118, e2016191118. \doi{10.1073/pnas.2016191118}
#'
#' \emph{Decomposition of proper scoring rules}
#'
#' Bröcker, J. (2009): `Reliability, sufficiency, and the decomposition of proper scores'. \emph{Quarterly Journal of the Royal Meteorological Society} 135, 1512-1519. \doi{10.1002/qj.456}
#'
#' \emph{Conditional decomposition of proper scoring rules}
#'
#' Allen, S., Ferro, C.A.T. and F. Kwasniok (2023): `A conditional decomposition of proper scores: quantifying the sources of information in a forecast'. \emph{Quarterly Journal of the Royal Meteorological Society}
#'
#'
#' @examples
#' p <- runif(100)
#' o <- rbinom(100, 1, 0.5)
#' p_disc <- sample(seq(0, 1, 0.1), 100, TRUE)
#'
#' # unconditional decomposition
#'
#' bs_decomp(o = o, p = p, bins = 10)
#' bs_decomp(o = o, p = p_disc)
#' bs_decomp(o = o, p = p, bins = 10, method = "classical")
#' bs_decomp(o = o, p = p, method = "isotonic")
#'
#' # conditional decomposition
#'
#' states = rep(c("G1", "G2"), each = 50)
#'
#' bs_decomp_cond(o = o, p = p, states = states, bins = 10)
#' bs_decomp_cond(o = o, p = p_disc, states = states)
#' bs_decomp_cond(o = o, p = p, states = states, bins = 10, method = "classical")
#' bs_decomp_cond(o = o, p = p, states = states, method = "isotonic")
#'
#' @name bs_decomp

#' @rdname bs_decomp
#' @export
bs_decomp <- function(o, p, bins = NULL, method = "isotonic") {
  check_input(o, p)
  if (method == "isotonic") {
    terms <- bs_decomp_iso(o, p)
  } else if (method %in% c("classical", "bias-corrected")) {

    if (is.null(bins)) {
      message("Taking bins to be the number of unique values of p")
      bins <- length(unique(p))
      if (bins/length(p) > 0.2) {
        warning("Number of bins is large relative to the number of observations.
                Consider specifying 'bins' manually.")
      }
    }

    p_breaks <- seq(0, 1, length.out = bins + 1)
    p_ind <- cut(p, breaks = p_breaks, include.lowest = TRUE, ordered_result = TRUE)
    p_disc <- (p_breaks[-(bins + 1)] + p_breaks[-1])/2
    p <- p_disc[p_ind]

    forecasts <- unique(p)
    obar <- mean(o)
    N <- length(o)
    n_k <- sapply(seq_along(forecasts), function(k) sum(p == forecasts[k]))
    obar_k <- sapply(seq_along(forecasts), function(k) mean(o[p == forecasts[k]]))

    unc <- obar*(1 - obar)
    res <- sum((n_k/N)*((obar_k - obar)^2))
    rel <- sum((n_k/N)*((forecasts - obar_k)^2))

    if(method == "bias-corrected"){
      if (any(n_k == 1)) {
        warning("There is a forecast value that only occurs once. The bias-corrected
                method cannot be used in this case. Returning estimates obtained without
                bias-correction.")
      } else {
        bias1 <- obar*(1 - obar)/(N - 1)
        bias2 <- sum(n_k*obar_k*(1 - obar_k)/(n_k - 1))/N
        unc <- unc + bias1
        res <- res + bias1 - bias2
        rel <- rel - bias2
      }
    }

    tot <- unc - res + rel
    terms <- c(unc, res, rel, tot)
    names(terms) <- c("UNC", "RES", "REL", "TOT")
  } else {
    stop("method not recognised. It must be one of 'classical', 'bias-corrected',
         and 'isotonic'")
  }
  return(terms)
}

#' @rdname bs_decomp
#' @export
bs_decomp_cond <- function(o, p, states, bins = NULL, method = "isotonic"){
  check_input(o, p, states)
  groups <- unique(states)
  n_j <- sapply(seq_along(groups), function(j) sum(states == groups[j]))
  N <- length(o)

  if(method == "isotonic"){

    terms_un <- bs_decomp_iso(o, p) # unconditional terms
    terms_cnd <- sapply(seq_along(groups), function(j)
      bs_decomp_iso(o[states == groups[j]], p[states == groups[j]])) # conditional terms

    unc_A <- sum((n_j/N)*terms_cnd['UNC', ])
    res_fA <- sum((n_j/N)*terms_cnd['RES', ])
    rel_fA <- sum((n_j/N)*terms_cnd['REL', ])

    res_A <- terms_un['UNC'] - unc_A
    res_Af <- rel_fA - terms_un['REL']

  }else if(method %in% c("classical", "bias-corrected")){

    if (is.null(bins)) {
      message("Taking bins to be the number of unique values of p")
      bins <- length(unique(p))
      if (bins/length(p) > 0.2) {
        warning("Number of bins is large relative to the number of observations.
                Consider specifying 'bins' manually.")
      }
    }

    p_breaks <- seq(0, 1, length.out = bins + 1)
    p_ind <- cut(p, breaks = p_breaks, include.lowest = TRUE, ordered_result = TRUE)
    p_disc <- (p_breaks[-(bins + 1)] + p_breaks[-1])/2
    p <- p_disc[p_ind]

    forecasts <- unique(p)
    obar <- mean(o)

    n_k <- sapply(seq_along(forecasts), function(k) sum(p == forecasts[k]))
    n_kj <- sapply(seq_along(forecasts), function(k)
      sapply(seq_along(groups), function(j) sum(p == forecasts[k] & states == groups[j])))

    obar_k <- sapply(seq_along(forecasts), function(k) mean(o[p == forecasts[k]]))
    obar_j <- sapply(seq_along(groups), function(j) mean(o[states == groups[j]]))
    obar_kj <- sapply(seq_along(forecasts), function(k)
      sapply(seq_along(groups), function(j) mean(o[p == forecasts[k] & states == groups[j]])))

    unc_A <- sum((n_j/N)*obar_j*(1 - obar_j))
    res_A <- sum((n_j/N)*((obar_j - obar)^2))
    res_fA <- sum((n_kj/N)*((obar_j - obar_kj)^2), na.rm = T)
    res_Af <- sum(t(n_kj/N)*((obar_k - t(obar_kj))^2), na.rm = T)
    rel_fA <- sum(t(n_kj/N)*((forecasts - t(obar_kj))^2), na.rm = T)
    tot <- unc_A - res_fA + rel_fA

    if(method == "bias-corrected"){
      if (any(n_k == 1)) {
        warning("There is a forecast value that only occurs once. The bias-corrected
                method cannot be used in this case. Returning estimates obtained without
                bias-correction.")
      } else if (any(n_j == 1)) {
        warning("There is a state that only occurs once. The bias-corrected
                method cannot be used in this case. Returning estimates obtained without
                bias-correction.")
      } else if (any(n_kj == 1)) {
        warning("There is a forecast-state combination that only occurs once. The bias-corrected
                method cannot be used in this case. Returning estimates obtained without
                bias-correction.")
      } else {
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
  } else {
    stop("method not recognised. It must be one of 'classical', 'bias-corrected',
         and 'isotonic'")
  }

  tot <- unc_A - res_fA + rel_fA
  terms <- c(unc_A, res_A, res_fA, res_Af, rel_fA, tot)
  names(terms) <- c("UNC_A", "RES_A", "RES_F|A", "RES_A|F", "REL_F|A", "TOT")

  return(terms)

}

# function to obtain Brier score decomposition terms from isotonic regression
bs_decomp_iso <- function(o, p, bins, states = NULL){
  diag_obj <- reliabilitydiag::reliabilitydiag(x = p, y = o)
  terms <- c(summary(diag_obj)$uncertainty,
             summary(diag_obj)$discrimination,
             summary(diag_obj)$miscalibration,
             summary(diag_obj)$mean_score)
  names(terms) <- c("UNC", "RES", "REL", "TOT")
  return(terms)
}

# function to check function inputs
check_input <- function(o, p, states = NULL) {
  if (!identical(length(o), length(p))) {
    stop("o and p must have the same length.")
  }
  if (any(is.na(o)) | any(is.na(p)) | any(is.na(states))) {
    stop("Input contains missing values.")
  }
  if (!all(o %in% c(0, 1))) {
    stop("o must only contain values that are either zero or one.")
  }
  if (!all(p >= 0 & p <= 1)) {
    stop("p must only contain values that are between zero and one.")
  }
  if (!is.null(states)) {
    if (!identical(length(o), length(states))) {
      stop("o, p and states must all have the same length.")
    }
  }
}




