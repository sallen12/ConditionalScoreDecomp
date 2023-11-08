#' Plots for score decompositions
#'
#' Plot decompositions of the proper scores into uncertainty, resolution, and
#' reliability components. A conditional decomposition can also be implemented, which
#' displays these terms conditional on some events, or states having occurred.
#'
#' @param terms_un vector containing the unconditional decomposition terms.
#' @param terms_cnd vector containing the conditional decomposition terms.
#' @param title optional title of the plot.
#' @param waterfall logical specifying whether a waterfall plot should be used.
#' @param dec_places integer specifying how many decimal places should be displayed.
#'
#' @details
#' These plot functions take decomposition terms as inputs, and return a plot
#' containing these values. Two plots are available, a bar plot (default) and a
#' waterfall plot (\code{waterfall = TRUE}).
#'
#' Decomposition terms should be inputted using \code{terms_un} and \code{terms_cnd}.
#' \code{terms_un} contains a vector with the uncertainty, resolution, reliability,
#' and total score. Note that this must be the ordering of the inputs.
#' \code{terms_unc} contains a vector with the conditional uncertainty, the resolution
#' of the states, the resolution of the forecasts given the states, the resolution
#' of the states given the forecasts, and the reliability of the forecasts given
#' the states.
#' These inputs are designed to correspond to the output from \code{\link{bs_decomp}}
#' and \code{\link{bs_decomp_cond}}.
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
#' scale <- 10000 # default plot parameters may not be suitable for all scales
#'
#' # unconditional decomposition
#' terms_un <- bs_decomp(o = o, p = p)*scale
#'
#' # conditional decomposition
#' states = rep(c("G1", "G2"), each = 500)
#' terms_cnd <- bs_decomp_cond(o = o, p = p, states = states)*scale
#'
#' # plot decomposition terms
#'
#' plot_decomp(terms_un)
#' plot_decomp(terms_un, dec_places = 0)
#' plot_decomp(terms_un, waterfall = TRUE)
#' plot_decomp(terms_un, waterfall = TRUE, dec_places = 0)
#'
#' plot_decomp(terms_cnd = terms_cnd)
#' plot_decomp(terms_cnd = terms_cnd, waterfall = TRUE)
#'
#' plot_decomp(terms_un, terms_cnd)
#' plot_decomp(terms_un, terms_cnd, waterfall = TRUE)
#'
#' @name plot_decomp

#' @rdname plot_decomp
#' @export
plot_decomp <- function(terms_un = NULL, terms_cnd = NULL, title = "", waterfall = FALSE, dec_places = 1){
  check_plot_input(terms_un, terms_cnd)
  if (waterfall) {
    plot_decomp_waterfall(terms_un, terms_cnd, title, dec_places)
  } else {
    plot_decomp_bar(terms_un, terms_cnd, title, dec_places)
  }
}

# function to plot decomposition terms in a waterfall plot
plot_decomp_waterfall <- function(terms_un = NULL, terms_cnd = NULL, title = "", dec_places = 1) {
  fmt <- paste("%.", dec_places, "f", sep = "")
  if (!is.null(terms_un)) {
    if (is.null(terms_cnd)) {
      ymin_id <- ymax_id <- c(1, 2, 3, 4)
    } else {
      ymin_id <- c(1, 3, 6, 8)
      ymax_id <- c(2, 5, 7, 8)
    }
    terms_un[2] <- -terms_un[2]
    df <- data.frame(begin = c(0, cumsum(terms_un[1:2]), 0),
                     end = terms_un + c(0, cumsum(terms_un[1:2]), 0),
                     id = (ymin_id + ymax_id)/2,
                     names = c("UNC_Y", "-RES_F", "REL_F", "Total"),
                     cols = as.factor(c(1, 2, 3, 4)))
    plot_un <- ggplot2::ggplot(df) +
      ggplot2::geom_rect(ggplot2::aes_string(ymin = "ymin_id - 0.45",
                                             ymax = "ymax_id + 0.45",
                                             xmin = "begin",
                                             xmax = "end",
                                             fill = "cols")) +
      ggplot2::geom_text(ggplot2::aes_string(x = "(begin + end)/2",
                                             y = "id",
                                             label = "sprintf(fmt, terms_un)"),
                         colour = c(rep("black", 3), "white")) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
      ggplot2::geom_linerange(ggplot2::aes_string(x = "end",
                                                  ymin = "ymax_id + 0.45",
                                                  ymax = "ymax_id + 1 - 0.45"),
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
    plot_cnd <- ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes_string(ymin = "id - 0.45",
                                                                             ymax = "id + 0.45",
                                                                             xmin = "begin",
                                                                             xmax = "end",
                                                                             fill = "cols")) +
      ggplot2::geom_text(ggplot2::aes_string(x = "(begin + end)/2",
                                             y = "id",
                                             label = "sprintf(fmt, terms_cnd)"),
                         colour = c(rep("black", 7), "white")) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
      ggplot2::geom_linerange(ggplot2::aes_string(x = "end",
                                                  ymin = "id + 0.45",
                                                  ymax = "id + 1 - 0.45"),
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
plot_decomp_bar <- function(terms_un = NULL, terms_cnd = NULL, title = "", dec_places = 1) {
  fmt <- paste("%.", dec_places, "f", sep = "")
  if (!is.null(terms_un)) {
    if (is.null(terms_cnd)) {
      ymin_id <- ymax_id <- c(1, 2, 3, 4)
    } else {
      ymin_id <- c(1, 3, 6, 8)
      ymax_id <- c(2, 5, 7, 8)
    }
    terms_un[2] <- -terms_un[2]
    df <- data.frame(begin = c(0, cumsum(terms_un[1:2]), 0),
                     end = terms_un + c(0, cumsum(terms_un[1:2]), 0),
                     id = (ymin_id + ymax_id)/2,
                     names = c(" UNC_Y", "-RES_F", " REL_F", " Total"),
                     cols = as.factor(c(1, 2, 3, 4)))
    lim_low <- min(terms_un) - 0.1*abs(max(terms_un))
    lim_upp <- max(terms_un) + 0.05*abs(max(terms_un))
    plot_un <- ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes_string(ymin = "ymin_id - 0.45",
                                                                            ymax = "ymax_id + 0.45",
                                                                            xmin = 0,
                                                                            xmax = "end - begin",
                                                                            fill = "cols")) +
      ggplot2::geom_text(ggplot2::aes_string(x = "lim_low + 0.075*abs(max(terms_un))",
                                      y = "id",
                                      label = "sprintf(fmt, terms_un)",
                                      fontface = 2),
                         colour = c(scales::hue_pal()(3), "black"), hjust = 1) +
      ggplot2::scale_fill_manual(values = c(scales::hue_pal()(3), "black")) +
      ggplot2::scale_x_continuous(name = NULL, limits = c(lim_low, lim_upp)) +
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
    terms_cnd <- c(terms_cnd[1], terms_cnd[2],
                   -terms_cnd[2], -terms_cnd[3], terms_cnd[4],
                   -terms_cnd[4], terms_cnd[5],
                   terms_cnd[6])
    if (is.null(terms_un)) {
      ylab_pos <- "left"
      title <- title
      lim_low <- min(terms_cnd) - 0.1*abs(max(terms_cnd))
      lim_upp <- max(terms_cnd) + 0.05*abs(max(terms_cnd))
      text_pos <- lim_low + 0.075*abs(max(terms_cnd))
    } else {
      ylab_pos <- "right"
      title <- ""
      lim_low <- min(terms_cnd) - 0.05*abs(max(terms_cnd))
      lim_upp <- max(terms_cnd) + 0.1*abs(max(terms_cnd))
      text_pos <- lim_upp - 0.025*abs(max(terms_cnd))
    }
  df <- data.frame(begin = c(0, cumsum(terms_cnd[1:6]), 0),
                   end = terms_cnd + c(0, cumsum(terms_cnd[1:6]), 0),
                   id = 1:8,
                   names = c("UNC_Y|A", "RES_A", "-RES_A", "-RES_F|A", "RES_A|F", "-RES_A|F", "REL_F|A", "Total"),
                   cols = as.factor(c(1, 1, 2, 2, 2, 3, 3, 4)))

  plot_cnd <- ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes_string(ymin = "id - 0.45",
                                                                           ymax = "id + 0.45",
                                                                           xmin = 0,
                                                                           xmax = "end - begin",
                                                                           fill = "cols")) +
    ggplot2::geom_text(ggplot2::aes_string(x = "text_pos",
                                           y = "id",
                                           label = "sprintf(fmt, terms_cnd)",
                                           fontface = 2),
                       colour = c(scales::hue_pal()(3), "black")[df$cols], hjust = 1) +
    ggplot2::scale_fill_manual(values = c(scales::hue_pal()(3), "black")) +
    ggplot2::scale_x_continuous(name = NULL, limits = c(lim_low, lim_upp)) +
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

# function to check plot inputs
check_plot_input <- function(terms_un, terms_cnd) {
  if (is.null(terms_un) && is.null(terms_cnd)) {
    stop("At least one of terms_un and terms_cnd must be specified.")
  }
  if (any(is.na(terms_un)) | any(is.na(terms_cnd))) {
    stop("Input contains missing values.")
  }
  if (!is.null(terms_un)) {
    if (!is.numeric(terms_un) | !is.vector(terms_un)) {
      stop("terms_un must be a vector of numeric values.")
    }
    if (!identical(length(terms_un), 4L)) {
      stop("terms_un must be of length 4, containing uncertainty, resolution,
         reliability and total score components.")
    }
  }
  if (!is.null(terms_cnd)) {
    if (!is.numeric(terms_cnd) | !is.vector(terms_cnd)) {
      stop("terms_cnd must be a vector of numeric values.")
    }
    if (!identical(length(terms_cnd), 6L)) {
      stop("terms_cnd must be of length 6.")
    }
  }
}
