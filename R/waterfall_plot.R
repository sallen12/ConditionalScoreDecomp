# function to create a waterfall plot
waterfall_plot <- function(p, o, states, title = "", scale = 10000){

  # plot parameters may need adjusting for scales different from 10000
  if(is.null(scale)){
    scale <- 1
  }

  # unconditional
  mth_cl <- bs_decomp(p, o)*scale
  mth_cl[2] <- -mth_cl[2]
  df <- data.frame(begin = c(0, cumsum(mth_cl[1:2]), 0), end = mth_cl + c(0, cumsum(mth_cl[1:2]), 0), id = c(1.5, 4, 6.5, 8),
                   names = c(" UNC_Y", "-RES_F", " REL_F", " Total"),
                   cols = as.factor(c(1, 2, 3, 4)))
  plot_cl <- ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes(ymin = c(1, 3, 6, 8) - 0.45, ymax = c(2, 5, 7, 8) + 0.45,
                                        xmin = 0, xmax = end - begin, fill = cols)) +
    ggplot2::geom_text(ggplot2::aes(x = 0.35*scale, y = id, label = sprintf("%.1f", mth_cl), fontface = 2),
              colour = c(scales::hue_pal()(3), "black"), hjust = 1) +
    scale_fill_manual(values = c(scales::hue_pal()(3), "black")) +
    scale_x_continuous(name = NULL, limits = c(-0.15, 0.35)*scale, breaks = seq(-0.1, 0.2, 0.1)*scale) +
    scale_y_continuous(name = NULL, breaks = df$id, position = "right",
                       labels = c(expression(paste("   ", UNC[Y])),
                                  expression(-RES["F"]),
                                  expression(paste("   ", REL["F"])),
                                  paste(" ", "Total"))) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
    ggplot2::ggtitle("") + theme_classic() +
    theme(legend.position = "none", axis.text = ggplot2::element_text(size = 11), axis.ticks.y = element_ggplot2::blank())

  # conditional
  mth <- bs_decomp_cond(p, o, states)*scale
  mth <- c(mth[1], mth[2], -mth[2], -mth[3], mth[4], -mth[4], mth[5], mth[6])
  df <- data.frame(begin = c(0, cumsum(mth[1:6]), 0), end = mth + c(0, cumsum(mth[1:6]), 0), id = 1:8,
                   names = c("UNC_Y|A", "RES_A", "-RES_A", "-RES_F|A", "RES_A|F", "-RES_A|F", "REL_F|A", "Total"),
                   cols = as.factor(c(1, 1, 2, 2, 2, 3, 3, 4)))
  plot_new <- ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes(ymin = id - 0.45, ymax = id + 0.45, xmin = 0, xmax = end - begin, fill = cols)) +
    ggplot2::geom_text(ggplot2::aes(x = -0.15*scale, y = id, label = sprintf("%.1f", mth), fontface = 2),
              colour = c(scales::hue_pal()(3), "black")[df$cols], hjust = 1) +
    scale_fill_manual(values = c(scales::hue_pal()(3), "black")) +
    scale_x_continuous(name = NULL, limits = c(-0.25, 0.25)*scale, breaks = seq(-0.1, 0.2, 0.1)*scale, expand = c(0, 0)) +
    scale_y_continuous(name = NULL, breaks = df$id,
                       labels = c(expression(paste("   ", UNC["Y|A"])),
                                  expression(paste("   ", RES[A])),
                                  expression(-RES[A]),
                                  expression(-RES["F|A"]),
                                  expression(paste("   ", RES["A|F"])),
                                  expression(-RES["A|F"]),
                                  expression(paste("   ", REL["F|A"])),
                                  paste("  ", "Total"))) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
    ggplot2::ggtitle(title) + theme_classic() +
    theme(legend.position = "none", axis.text = ggplot2::element_text(size = 11),
          axis.ticks.y = element_ggplot2::blank(), axis.text.y = ggplot2::element_text(hjust = 0))

  gridExtra::grid.arrange(plot_new, plot_cl, nrow = 1)
}


# function to create an alternative waterfall plot
waterfall_plot_alt <- function(p, o, states, title = "", scale = 10000){

  # plot parameters may need adjusting for scales different from 10000
  if(is.null(scale)){
    scale <- 1
  }

  # unconditional
  mth_cl <- bs_decomp(p, o)*scale
  mth_cl[2] <- -mth_cl[2]
  df <- data.frame(begin = c(0, cumsum(mth_cl[1:2]), 0), end = mth_cl + c(0, cumsum(mth_cl[1:2]), 0),
                   id = c(1.5, 4, 6.5, 8),
                   names = c("UNC_Y", "-RES_F", "REL_F", "Total"),
                   cols = as.factor(c(1, 2, 3, 4)))
  plot_cl <- ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes(ymin = c(1, 3, 6, 8) - 0.45, ymax = c(2, 5, 7, 8) + 0.45,
                                        xmin = begin, xmax = end, fill = cols)) +
    ggplot2::geom_text(ggplot2::aes(x = (begin + end)/2, y = id, label = sprintf("%.1f", mth_cl)), colour = c(rep("black", 3), "white")) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
    ggplot2::geom_linerange(ggplot2::aes(x = end, ymin = c(2, 5, 7, 8) + 0.45, ymax = c(2, 5, 7, 8) + 1 - 0.45),
                   size = 1, colour = c(rep("black", 3), "white")) +
    scale_fill_manual(values = c(scales::hue_pal()(3), "black")) +
    scale_x_continuous(name = NULL, limits = c(0, 2600), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, breaks = df$id, position = "right",
                       labels = c(expression(UNC[Y]), expression(-RES["F"]), expression(REL["F"]), "Total")) +
    ggplot2::ggtitle("") + theme_classic() +
    theme(legend.position = "none", axis.text = ggplot2::element_text(size = 11))

  # conditional
  mth <- bs_decomp_cond(p, o, states)*scale
  mth <- c(mth[1], mth[2], -mth[2], -mth[3], mth[4], -mth[4], mth[5], mth[6])
  df <- data.frame(begin = c(0, cumsum(mth[1:6]), 0), end = mth + c(0, cumsum(mth[1:6]), 0), id = 1:8,
                   names = c("UNC_Y|A", "RES_A", "-RES_A", "-RES_F|A", "RES_A|F", "-RES_A|F", "REL_F|A", "Total"),
                   cols = as.factor(c(1, 1, 2, 2, 2, 3, 3, 4)))
  plot_new <- ggplot2::ggplot(df) + ggplot2::geom_rect(ggplot2::aes(ymin = id - 0.45, ymax = id + 0.45, xmin = begin, xmax = end, fill = cols)) +
    ggplot2::geom_text(ggplot2::aes(x = (begin + end)/2, y = id, label = sprintf("%.1f", mth)), colour = c(rep("black", 7), "white")) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0)) +
    ggplot2::geom_linerange(ggplot2::aes(x = end, ymin = id + 0.45, ymax = id + 1 - 0.45),
                   size = 1, colour = c(rep("black", 7), "white")) +
    scale_fill_manual(values = c(scales::hue_pal()(3), "black")) +
    scale_x_continuous(name = NULL, limits = c(0, 2600), expand = c(0, 0)) +
    scale_y_continuous(name = NULL, breaks = df$id,
                       labels = c(expression(UNC["Y|A"]), expression(RES[A]), expression(-RES[A]),
                                  expression(-RES["F|A"]), expression(RES["A|F"]), expression(-RES["A|F"]),
                                  expression(REL["F|A"]), "Total")) +
    ggplot2::ggtitle(title) + theme_classic() +
    theme(legend.position = "none", axis.text = ggplot2::element_text(size = 11))

  gridExtra::grid.arrange(plot_new, plot_cl, nrow = 1)
}
