#' Plots for score decompositions
#'
#' Plot decompositions of the Brier score into uncertainty, resolution, and
#' reliability components. A conditional decomposition can also be plot, which
#' displays these terms conditional on some events, or states having occurred.
#'
#' @param terms_un vector containing the unconditional decomposition terms.
#' @param terms_cnd vector containing the conditional decomposition terms.
#' @param title optional title of the plot.
#' @param scale amount by which to scale the decomposition terms.
#' @param waterfall logical specifying whether a waterfall plot should be used.
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
#' for small sample sizes; \code{method = "isotonic"} performs the
#' decomposition based on isotonic regression proposed by Dimitriadies et al. (2021).
#' See references below.
#'
#' Calculation of the unconditional Brier score decomposition (\code{bs_decomp})
#' leverages the \pkg{SpecsVerification} and \pkg{reliabilitydiag} packages.
#'
#' @author Sam Allen
#'
#' @references
#'
#' Allen, S., Ferro, C.A.T. and F. Kwasniok (2023): `A conditional decomposition of proper scores: quantifying the sources of information in a forecast'. \emph{Quarterly Journal of the Royal Meteorological Society}
#'
#'
#' @examples
#' o <- rbinom(1000, 1, 0.5)
#' p <- sample(seq(0, 1, 0.1), 1000, TRUE)
#'
#' # unconditional decomposition
#' terms_un <- bs_decomp(o = o, p = p)
#'
#' # conditional decomposition
#' states = rep(c("G1", "G2"), each = 500)
#' terms_cnd <- bs_decomp_cond(o = o, p = p, states = states)
#'
#' # plot decomposition terms
#'
#' plot_decomp(terms_un)
#' plot_decomp(terms_un, waterfall = TRUE)
#'
#' plot_decomp(terms_cnd = terms_cnd)
#' plot_decomp(terms_cnd = terms_cnd, waterfall = TRUE)
#'
#' plot_decomp(terms_un, terms_cnd)
#' plot_decomp(terms_un, terms_cnd, waterfall = TRUE)
#'
#'
#' @name plot_decomp

#' @rdname plot_decomp
#' @export
plot_decomp <- function(terms_un = NULL, terms_cnd = NULL, title = "", waterfall = FALSE){
  if (waterfall) {
    plot_decomp_waterfall(terms_un, terms_cnd, title)
  } else {
    plot_decomp_bar(terms_un, terms_cnd, title = title)
  }
}

# function to plot decomposition terms in a waterfall plot
plot_decomp_waterfall <- function(terms_un = NULL, terms_cnd = NULL, title = "") {

  if (!is.null(terms_un)) {
    terms_un[2] <- -terms_un[2]
    df <- data.frame(begin = c(0, cumsum(terms_un[1:2]), 0),
                     end = terms_un + c(0, cumsum(terms_un[1:2]), 0),
                     id = c(1.5, 4, 6.5, 8),
                     names = c("UNC_Y", "-RES_F", "REL_F", "Total"),
                     cols = as.factor(c(1, 2, 3, 4)))
    plot_un <- ggplot2::ggplot(df) +
      ggplot2::geom_rect(ggplot2::aes(ymin = c(1, 3, 6, 8) - 0.45,
                                      ymax = c(2, 5, 7, 8) + 0.45,
                                      xmin = begin,
                                      xmax = end,
                                      fill = cols)) +
      ggplot2::geom_text(ggplot2::aes(x = (begin + end)/2,
                                      y = id,
                                      label = sprintf("%.1f", terms_un)),
                         colour = c(rep("black", 3), "white")) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
      ggplot2::geom_linerange(ggplot2::aes(x = end, ymin = c(2, 5, 7, 8) + 0.45, ymax = c(2, 5, 7, 8) + 1 - 0.45),
                              linewidth = 1, colour = c(rep("black", 3), "white")) +
      ggplot2::scale_fill_manual(values = c(scales::hue_pal()(3), "black")) +
      ggplot2::scale_x_continuous(name = NULL, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(name = NULL, breaks = df$id,
                                  labels = c(expression(UNC[Y]), expression(-RES["F"]), expression(REL["F"]), "Total")) +
      ggplot2::ggtitle(title) + ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none", axis.text = ggplot2::element_text(size = 11))
  }

  if (!is.null(terms_cnd)) {
    if (is.null(terms_un)) {
      ylab_pos <- "left"
      title <- title
    } else {
      ylab_pos <- "right"
      title <- ""
    }
    terms_cnd <- c(terms_cnd[1], terms_cnd[2],
                   -terms_cnd[2], -terms_cnd[3], terms_cnd[4],
                   -terms_cnd[4], terms_cnd[5],
                   terms_cnd[6])
    df <- data.frame(begin = c(0, cumsum(terms_cnd[1:6]), 0),
                     end = terms_cnd + c(0, cumsum(terms_cnd[1:6]), 0),
                     id = 1:8,
                     names = c("UNC_Y|A", "RES_A", "-RES_A", "-RES_F|A", "RES_A|F", "-RES_A|F", "REL_F|A", "Total"),
                     cols = as.factor(c(1, 1, 2, 2, 2, 3, 3, 4)))
    plot_cnd <- ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes(ymin = id - 0.45,
                                                                      ymax = id + 0.45,
                                                                      xmin = begin,
                                                                      xmax = end,
                                                                      fill = cols)) +
      ggplot2::geom_text(ggplot2::aes(x = (begin + end)/2,
                                      y = id,
                                      label = sprintf("%.1f", terms_cnd)),
                         colour = c(rep("black", 7), "white")) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
      ggplot2::geom_linerange(ggplot2::aes(x = end,
                                           ymin = id + 0.45,
                                           ymax = id + 1 - 0.45),
                              linewidth = 1,
                              colour = c(rep("black", 7), "white")) +
      ggplot2::scale_fill_manual(values = c(scales::hue_pal()(3), "black")) +
      ggplot2::scale_x_continuous(name = NULL, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(name = NULL, breaks = df$id, position = ylab_pos,
                                  labels = c(expression(UNC["Y|A"]), expression(RES[A]), expression(-RES[A]),
                                             expression(-RES["F|A"]), expression(RES["A|F"]), expression(-RES["A|F"]),
                                             expression(REL["F|A"]), "Total")) +
      ggplot2::ggtitle(title) + ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none", axis.text = ggplot2::element_text(size = 11))
  }

  if (is.null(terms_cnd)) {
    plot_un
  } else if (is.null(terms_un)) {
    plot_cnd
  } else {
    gridExtra::grid.arrange(plot_un, plot_cnd, nrow = 1)
  }
}

# function to plot decomposition terms in a bar plot
plot_decomp_bar <- function(terms_un = NULL, terms_cnd = NULL, title = "", scale = 10000) {

  if (!is.null(terms_un)) {
    terms_un[2] <- -terms_un[2]
    df <- data.frame(begin = c(0, cumsum(terms_un[1:2]), 0),
                     end = terms_un + c(0, cumsum(terms_un[1:2]), 0),
                     id = c(1.5, 4, 6.5, 8),
                     names = c(" UNC_Y", "-RES_F", " REL_F", " Total"),
                     cols = as.factor(c(1, 2, 3, 4)))
    plot_un <- ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes(ymin = c(1, 3, 6, 8) - 0.45,
                                                                     ymax = c(2, 5, 7, 8) + 0.45,
                                                                     xmin = 0,
                                                                     xmax = end - begin,
                                                                     fill = cols)) +
      ggplot2::geom_text(ggplot2::aes(x = min(terms_un),
                                      y = id,
                                      label = sprintf("%.1f", terms_un),
                                      fontface = 2),
                         colour = c(scales::hue_pal()(3), "black"), hjust = 1) +
      ggplot2::scale_fill_manual(values = c(scales::hue_pal()(3), "black")) +
      ggplot2::scale_x_continuous(name = NULL) +
      ggplot2::scale_y_continuous(name = NULL, breaks = df$id,
                                  labels = c(expression(paste("   ", UNC[Y])),
                                             expression(-RES["F"]),
                                             expression(paste("   ", REL["F"])),
                                             paste(" ", "Total"))) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
      ggplot2::ggtitle(title) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none",
                     axis.text = ggplot2::element_text(size = 11),
                     axis.ticks.y = ggplot2::element_blank())
  }

  if (!is.null(terms_cnd)) {
    if (is.null(terms_un)) {
      ylab_pos <- "left"
      title <- title
    } else {
      ylab_pos <- "right"
      title <- ""
    }
    terms_cnd <- c(terms_cnd[1], terms_cnd[2],
                   -terms_cnd[2], -terms_cnd[3], terms_cnd[4],
                   -terms_cnd[4], terms_cnd[5],
                   terms_cnd[6])
  df <- data.frame(begin = c(0, cumsum(terms_cnd[1:6]), 0),
                   end = terms_cnd + c(0, cumsum(terms_cnd[1:6]), 0),
                   id = 1:8,
                   names = c("UNC_Y|A", "RES_A", "-RES_A", "-RES_F|A", "RES_A|F", "-RES_A|F", "REL_F|A", "Total"),
                   cols = as.factor(c(1, 1, 2, 2, 2, 3, 3, 4)))
  plot_cnd <- ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes(ymin = id - 0.45,
                                                                    ymax = id + 0.45,
                                                                    xmin = 0,
                                                                    xmax = end - begin,
                                                                    fill = cols)) +
    ggplot2::geom_text(ggplot2::aes(x = max(terms_cnd),
                                    y = id,
                                    label = sprintf("%.1f", terms_cnd),
                                    fontface = 2),
                       colour = c(scales::hue_pal()(3), "black")[df$cols], hjust = 1) +
    ggplot2::scale_fill_manual(values = c(scales::hue_pal()(3), "black")) +
    ggplot2::scale_x_continuous(name = NULL, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = NULL, breaks = df$id, position = ylab_pos,
                                labels = c(expression(paste("   ", UNC["Y|A"])),
                                           expression(paste("   ", RES[A])),
                                           expression(-RES[A]),
                                           expression(-RES["F|A"]),
                                           expression(paste("   ", RES["A|F"])),
                                           expression(-RES["A|F"]),
                                           expression(paste("   ", REL["F|A"])),
                                           paste("  ", "Total"))) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
    ggplot2::ggtitle(title) + ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text = ggplot2::element_text(size = 11),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(hjust = 0))
  }

  if (is.null(terms_cnd)) {
    plot_un
  } else if (is.null(terms_un)) {
    plot_cnd
  } else {
    gridExtra::grid.arrange(plot_un, plot_cnd, nrow = 1)
  }
}
