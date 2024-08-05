options(width = 300)

library(tidyverse)
library(latex2exp)
library(RColorBrewer)

# APPENDIX ---------------------------------------------------------------------
## SIMPLIFIED METHODS ----------------------------------------------------------
### Block randomized response model  -------------------------------------------
### Plot of the finite sample correction for the optimized betas and the rr betas
load_data <- function(exp.num) {
  idir <- sprintf("results/exp%d/simplified_methods", exp.num)
  ifile.list <- list.files(idir)
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  
  
  #summary <- results
  summary <- results %>%
    pivot_longer("n_cal", names_to = "Key", values_to = "Value")# %>%
  #group_by(data, num_var, K, signal, model_name, contamination, epsilon, estimate, n_train, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
  #summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
  
  return(results)
}

exp.num <- 1
summary <- load_data(exp.num)

# prova <- summary  %>%
#   pivot_longer("n_cal", names_to = "Key", values_to = "Value")

method.values = c("RR", "CVX")
method.labels = c("RR", "Opt")
cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
color.scale <- cbPalette[c(3,7)]
shape.scale <- c(2,0)
linetype.scale <- c(1,1)

#' Plot 1: Finite sample correction as a function of the calibration set size,
#'         for different number of classes
make_figure_1 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="Klab", epsilon==plot.epsilon, model_name=="B") %>%
    filter(K %in% c(4, 8, 16)) %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    mutate(K = factor(K, levels = label.values, labels = label.labels))
  
  pp <- ggplot(df, aes(x = n_cal, y = values, color = Method, shape = Method, linetype = Method)) +
    geom_point() +
    geom_line() +
    facet_wrap(~ K, scales = "fixed") +
    #scale_x_log10(breaks = unique(df$n_cal), labels = unique(df$n_cal)) +
    scale_x_continuous(trans = 'log10', breaks = unique(df$n_cal), labels = unique(df$n_cal)) +
    scale_color_manual(values = color.scale) +
    scale_shape_manual(values = shape.scale) +
    scale_linetype_manual(values = linetype.scale) +
    xlab("Number of calibration samples") +
    ylab(expression(delta^"marg" * "(n)")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/delta_marg_const_B_Klab_eps%f.pdf", plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(4, 8, 16)
label.labels = c("K=4", "K=8", "K=16")

make_figure_1(plot.epsilon=0.1, save_plots=TRUE, reload=FALSE)

#' Plot 2: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size
make_figure_2 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="nlab", epsilon==plot.epsilon, model_name=="B") %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    mutate(n_lab = factor(n_cal, levels = label.values, labels = label.labels))
  
  pp <- ggplot(df, aes(x = K, y = values, color = Method, shape = Method, linetype = Method)) +
    geom_point() +
    geom_line() +
    facet_wrap(~ n_lab, scales = "fixed") +
    scale_color_manual(values = color.scale) +
    scale_shape_manual(values = shape.scale) +
    scale_linetype_manual(values = linetype.scale) +
    xlab("Number of classes") +
    ylab(expression(delta^"marg" * "(n)")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/delta_marg_const_B_nlab_eps%f.pdf", plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(50, 100, 500)
label.labels = c("n=50", "n=100", "n=500")

make_figure_2(plot.epsilon=0.1, save_plots=TRUE, reload=FALSE)

### Two-level randomized response model ----------------------------------------
### Plot of the finite sample correction for the optimized betas and the rr betas
load_data <- function(exp.num) {
  idir <- sprintf("results/exp%d/simplified_methods", exp.num)
  ifile.list <- list.files(idir)
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  
  
  #summary <- results
  summary <- results %>%
    pivot_longer("n_cal", names_to = "Key", values_to = "Value")# %>%
    #group_by(data, num_var, K, signal, model_name, contamination, epsilon, estimate, n_train, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    #summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
  
  return(results)
}

exp.num <- 1
summary <- load_data(exp.num)

# prova <- summary  %>%
#   pivot_longer("n_cal", names_to = "Key", values_to = "Value")

method.values = c("RR", "CVX")
method.labels = c("RR", "Opt")
cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
color.scale <- cbPalette[c(3,7)]
shape.scale <- c(2,0)
linetype.scale <- c(1,1)

#' Plot 1: Finite sample correction as a function of the calibration set size,
#'         for different number of classes
make_figure_1 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="Klab", epsilon==plot.epsilon, nu==plot.nu, model_name=="BRR") %>%
    filter(K %in% c(4, 8, 16)) %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    mutate(K = factor(K, levels = label.values, labels = label.labels))
  
  pp <- ggplot(df, aes(x = n_cal, y = values, color = Method, shape = Method, linetype = Method)) +
    geom_point() +
    geom_line() +
    facet_wrap(~ K, scales = "fixed") +
    #scale_x_log10(breaks = unique(df$n_cal), labels = unique(df$n_cal)) +
    scale_x_continuous(trans = 'log10', breaks = unique(df$n_cal), labels = unique(df$n_cal)) +
    scale_color_manual(values = color.scale) +
    scale_shape_manual(values = shape.scale) +
    scale_linetype_manual(values = linetype.scale) +
    xlab("Number of calibration samples") +
    ylab(expression(delta^"marg" * "(n)")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/delta_marg_const_BRR_Klab_eps%f_nu%f.pdf", plot.epsilon, plot.nu)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(4, 8, 16)
label.labels = c("K=4", "K=8", "K=16")

make_figure_1(plot.epsilon=0.5, plot.nu=0.5, save_plots=TRUE, reload=FALSE)

#' Plot 2: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size
make_figure_2 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="nlab", epsilon==plot.epsilon, nu==plot.nu, model_name=="BRR") %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    mutate(n_lab = factor(n_cal, levels = label.values, labels = label.labels))
  
  pp <- ggplot(df, aes(x = K, y = values, color = Method, shape = Method, linetype = Method)) +
    geom_point() +
    geom_line() +
    facet_wrap(~ n_lab, scales = "fixed") +
    scale_color_manual(values = color.scale) +
    scale_shape_manual(values = shape.scale) +
    scale_linetype_manual(values = linetype.scale) +
    xlab("Number of classes") +
    ylab(expression(delta^"marg" * "(n)")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/delta_marg_const_BRR_nlab_eps%f_nu%f.pdf", plot.epsilon, plot.nu)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(50, 100, 500)
label.labels = c("n=50", "n=100", "n=500")

make_figure_2(plot.epsilon=0.5, plot.nu=0.5, save_plots=TRUE, reload=FALSE)


#' Plot 3: Finite sample correction as a function of the calibration set size,
#'         for different number of classes and for different combinations of epsilon and nu
#'

make_figure_3 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon.vals=c(0.1, 0.5), plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="Klab", epsilon %in% plot.epsilon.vals, nu==plot.nu, model_name=="BRR") %>%
    filter(K %in% c(4, 8, 16)) %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    mutate(K = factor(K, levels = label.values, labels = label.labels)) %>%
    mutate(epsilon = factor(epsilon, labels = paste0("epsilon=", plot.epsilon.vals)))
  
  y_min <- min(df$values, na.rm = TRUE)
  y_max <- max(df$values, na.rm = TRUE)
  
  pp <- ggplot(df, aes(x = n_cal, y = values, color = Method, shape = Method, linetype = Method)) +
    geom_point() +
    geom_line() +
    facet_grid(epsilon ~ K, scales = "fixed") +
    #scale_x_log10(breaks = unique(df$n_cal), labels = unique(df$n_cal)) +
    scale_x_continuous(trans = 'log10', breaks = unique(df$n_cal), labels = unique(df$n_cal)) +
    scale_color_manual(values = color.scale) +
    scale_shape_manual(values = shape.scale) +
    scale_linetype_manual(values = linetype.scale) +
    xlab("Number of calibration samples") +
    ylab(expression(delta^"marg" * "(n)")) +
    ylim(y_min, y_max) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/delta_marg_const_BRR_Klab_nu%f.pdf", plot.nu)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(4, 8, 16)
label.labels = c("K=4", "K=8", "K=16")

make_figure_3(plot.epsilon.vals=c(0.1,0.5), plot.nu=0.2, save_plots=TRUE, reload=FALSE)


        
#' Plot 4: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size and for different combinations of epsilon and nu
#'         
make_figure_4 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon.vals=c(0.1,0.5), plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="nlab", epsilon %in% plot.epsilon.vals, nu==plot.nu, model_name=="BRR") %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    mutate(n_lab = factor(n_cal, levels = label.values, labels = label.labels)) %>%
    mutate(epsilon = factor(epsilon, labels = paste0("epsilon=", plot.epsilon.vals)))
  
  y_min <- min(df$values, na.rm = TRUE)
  y_max <- max(df$values, na.rm = TRUE)
  
  pp <- ggplot(df, aes(x = K, y = values, color = Method, shape = Method, linetype = Method)) +
    geom_point() +
    geom_line() +
    facet_grid(epsilon ~ n_lab, scales = "fixed") +
    scale_color_manual(values = color.scale) +
    scale_shape_manual(values = shape.scale) +
    scale_linetype_manual(values = linetype.scale) +
    xlab("Number of classes") +
    ylab(expression(delta^"marg" * "(n)")) +
    ylim(y_min, y_max) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/delta_marg_const_BRR_nlab_nu%f.pdf", plot.nu)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(50, 100, 500)
label.labels = c("n=50", "n=100", "n=500")

make_figure_4(plot.epsilon.vals=c(0.1,0.5), plot.nu=0.2, save_plots=TRUE, reload=FALSE)

#' Plot 5: Finite sample correction as a function of the calibration set size,
#'         for different number of classes and for different combinations of epsilon and nu
#'

make_figure_5 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8)) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="Klab", epsilon == plot.epsilon, nu %in% plot.nu.vals, model_name=="BRR") %>%
    filter(K %in% c(4, 8, 16)) %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    mutate(K = factor(K, levels = label.values, labels = label.labels)) %>%
    mutate(nu = factor(nu, labels = paste0("nu=", plot.nu.vals)))

  y_min <- min(df$values, na.rm = TRUE)
  y_max <- max(df$values, na.rm = TRUE)
  
  pp <- ggplot(df, aes(x = n_cal, y = values, color = Method, shape = Method, linetype = Method)) +
    geom_point() +
    geom_line() +
    facet_grid(nu ~ K, scales = "fixed") +
    #scale_x_log10(breaks = unique(df$n_cal), labels = unique(df$n_cal)) +
    scale_x_continuous(trans = 'log10', breaks = unique(df$n_cal), labels = unique(df$n_cal)) +
    scale_color_manual(values = color.scale) +
    scale_shape_manual(values = shape.scale) +
    scale_linetype_manual(values = linetype.scale) +
    xlab("Number of calibration samples") +
    ylab(expression(delta^"marg" * "(n)")) +
    ylim(y_min, y_max) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/delta_marg_const_BRR_Klab_eps%f.pdf", plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(4, 8, 16)
label.labels = c("K=4", "K=8", "K=16")

make_figure_5(plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8), save_plots=TRUE, reload=FALSE)



#' Plot 6: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size and for different combinations of epsilon and nu
#'         
make_figure_6 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8)) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="nlab", epsilon == plot.epsilon, nu %in% plot.nu.vals, model_name=="BRR") %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    mutate(n_lab = factor(n_cal, levels = label.values, labels = label.labels)) %>%
    mutate(nu = factor(nu, labels = paste0("nu=", plot.nu.vals)))
  
  y_min <- min(df$values, na.rm = TRUE)
  y_max <- max(df$values, na.rm = TRUE)
  
  pp <- ggplot(df, aes(x = K, y = values, color = Method, shape = Method, linetype = Method)) +
    geom_point() +
    geom_line() +
    facet_grid(nu ~ n_lab, scales = "fixed") +
    scale_color_manual(values = color.scale) +
    scale_shape_manual(values = shape.scale) +
    scale_linetype_manual(values = linetype.scale) +
    xlab("Number of classes") +
    ylab(expression(delta^"marg" * "(n)")) +
    ylim(y_min, y_max) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/delta_marg_const_BRR_nlab_eps%f.pdf", plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(50, 100, 500)
label.labels = c("n=50", "n=100", "n=500")

make_figure_6(plot.epsilon=0.1, plot.nu.vals=c(0.2, 0.8), save_plots=TRUE, reload=FALSE)


# MAIN -------------------------------------------------------------------------


