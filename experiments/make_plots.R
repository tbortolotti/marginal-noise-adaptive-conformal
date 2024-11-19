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
  ifile.list <- list.files(idir, recursive=FALSE)
  ifile.list <- ifile.list[!grepl("comparison", ifile.list)]
  
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
make_figure_A1 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, model_name="B",
                           method.values, method.labels, label.values, label.labels) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="Klab", epsilon==plot.epsilon, model_name==model_name, K %in% label.values) %>%
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
    plot.file <- sprintf("figures/delta_marg_const_%s_Klab_eps%f.pdf", model_name, plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(4, 8, 16)
label.labels = c("K=4", "K=8", "K=16")
model.name = "B"
plot.epsilon = 0.1

make_figure_A1(plot.epsilon=plot.epsilon, model_name=model_name, method.values=method.values, method.labels=method.labels,
               label.values=label.values, label.labels=label.labels, save_plots=TRUE, reload=FALSE)

#' Plot 2: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size
make_figure_A2 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, model_name="B") {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="nlab", epsilon==plot.epsilon, model_name==model_name, n_cal %in% label.values) %>%
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
    plot.file <- sprintf("figures/delta_marg_const_%s_nlab_eps%f.pdf", model_name, plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(100, 500, 1000)
label.labels = c("n=100", "n=500", "n=1000")
plot.epsilon = 0.1
model_name = "B"

make_figure_A2(plot.epsilon=plot.epsilon, model_name=model_name, save_plots=TRUE, reload=FALSE)

### Two-level randomized response model ----------------------------------------
### Plot of the finite sample correction for the optimized betas and the rr betas
load_data <- function(exp.num) {
  idir <- sprintf("results/exp%d/simplified_methods", exp.num)
  ifile.list <- list.files(idir, recursive=FALSE)
  ifile.list <- ifile.list[!grepl("comparison", ifile.list)]
  
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

#' Plot 3: Finite sample correction as a function of the calibration set size,
#'         for different number of classes
make_figure_A3 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu=0.2, model_name="B") {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="Klab", epsilon==plot.epsilon, nu==plot.nu, model_name==model_name) %>%
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
    plot.file <- sprintf("figures/delta_marg_const_%s_Klab_eps%f_nu%f.pdf", model_name, plot.epsilon, plot.nu)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(4, 8, 16)
label.labels = c("K=4", "K=8", "K=16")
plot.epsilon = 0.1
plot.nu = 0.2
model_name = "BRR"

make_figure_A3(plot.epsilon=plot.epsilon, plot.nu=plot.nu, model_name=model_name, save_plots=TRUE, reload=FALSE)

#' Plot 4: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size
make_figure_A4 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu=0.2, model_name="B") {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="nlab", epsilon==plot.epsilon, nu==plot.nu, model_name==model_name,
           n_cal %in% label.values) %>%
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
    plot.file <- sprintf("figures/delta_marg_const_%s_nlab_eps%f_nu%f.pdf", model_name, plot.epsilon, plot.nu)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(100, 500, 1000)
label.labels = c("n=100", "n=500", "n=1000")
plot.epsilon = 0.1
plot.nu = 0.2
model_name = "BRR"

make_figure_A4(plot.epsilon=plot.epsilon, plot.nu=plot.nu, model_name=model_name, save_plots=TRUE, reload=FALSE)


#' Plot 5: Finite sample correction as a function of the calibration set size,
#'         for different number of classes and for different combinations of epsilon and nu
#'

make_figure_A5 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8)) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="Klab", epsilon == plot.epsilon, nu %in% plot.nu.vals, model_name=="BRR") %>%
    filter(K %in% label.values) %>%
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

make_figure_A5(plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8), save_plots=TRUE, reload=FALSE)



#' Plot 6: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size and for different combinations of epsilon and nu
#'         
make_figure_A6 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8)) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="nlab", epsilon == plot.epsilon, nu %in% plot.nu.vals, model_name=="BRR",
           n_cal %in% label.values) %>%
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

label.values = c(100, 500, 1000)
label.labels = c("n=100", "n=500", "n=1000")

make_figure_A6(plot.epsilon=0.1, plot.nu.vals=c(0.2, 0.8), save_plots=TRUE, reload=FALSE)

#' Plot 7: Finite sample correction as a function of nu, for different
#'         values of the calibration set size and for different number of classes
#'         
make_figure_A7 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.K.vals=c(4,16), model_name) {
  if(reload) {
    summary <- load_data(1)
  }
  df <- summary %>%
    filter(plot=="nu_var", epsilon == plot.epsilon, K %in% plot.K.vals, model_name==model_name,
           n_cal %in% label.values) %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    mutate(n_lab = factor(n_cal, levels = label.values, labels = label.labels)) %>%
    mutate(K = factor(K, levels=plot.K.vals, labels = c("K=4","K=16")))
  
  y_min <- min(df$values, na.rm = TRUE)
  y_max <- max(df$values, na.rm = TRUE)
  
  pp <- ggplot(df, aes(x = nu, y = values, color = Method, shape = Method, linetype = Method)) +
    geom_point() +
    geom_line() +
    facet_grid(K ~ n_lab, scales = "fixed") +
    scale_color_manual(values = color.scale) +
    scale_shape_manual(values = shape.scale) +
    scale_linetype_manual(values = linetype.scale) +
    xlab(expression(nu)) +
    ylab(expression(delta^"marg" * "(n)")) +
    ylim(y_min, y_max) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/delta_marg_const_BRR_nu_var_eps%f.pdf", plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(100, 500, 1000)
label.labels = c("n=100", "n=500", "n=1000")
#plot.epsilon = 0.1
plot.epsilon = 0.5
plot.K.vals = c(4,16)
model_name = "BRR"

make_figure_A7(plot.epsilon=plot.epsilon, plot.K.vals=plot.K.vals, model_name=model_name,
               save_plots=TRUE, reload=FALSE)



### Comparison with the old finite sample correction ----------------------------------------
### Plot of the finite sample correction for the optimized betas, the rr betas and the old
### finite sample correction
load_data <- function(exp.num) {
  idir <- sprintf("results/exp%d/simplified_methods/comparison", exp.num)
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

method.values = c("RR", "CVX", "OLD")
method.labels = c("RR", "Opt", "Old")
cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
color.scale <- cbPalette[c(3,7,8)]
shape.scale <- c(2,0,1)
linetype.scale <- c(1,1,1)

#' Plot 9: Finite sample correction as a function of the calibration set size,
#'         for different number of classes
make_figure_A9 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="Klab", epsilon==plot.epsilon, nu==plot.nu, model_name=="block") %>%
    filter(K %in% label.values) %>%
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
    plot.file <- sprintf("figures/comparison_delta_marg_const_block_Klab_eps%f.pdf", plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(4, 8, 16)
label.labels = c("K=4", "K=8", "K=16")

make_figure_A9(plot.epsilon=0.1, plot.nu="none", save_plots=TRUE, reload=FALSE)

#' Plot 10: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size
make_figure_A10 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="nlab", epsilon==plot.epsilon, nu==plot.nu, model_name=="block",
           n_cal %in% label.values) %>%
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
    plot.file <- sprintf("figures/comparison_delta_marg_const_block_nlab_eps%f.pdf", plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(500, 1000, 2000)
label.labels = c("n=500", "n=1000", "n=2000")

make_figure_A10(plot.epsilon=0.1, plot.nu="none", save_plots=TRUE, reload=FALSE)


#' Plot 11: Finite sample correction as a function of the calibration set size,
#'         for different number of classes
make_figure_A11 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="Klab", epsilon==plot.epsilon, nu==plot.nu, model_name=="BRR") %>%
    filter(K %in% label.values) %>%
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
    plot.file <- sprintf("figures/comparison_delta_marg_const_BRR_Klab_eps%f_nu%f.pdf", plot.epsilon, plot.nu)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(4, 8, 16)
label.labels = c("K=4", "K=8", "K=16")

make_figure_A11(plot.epsilon=0.1, plot.nu=0.2, save_plots=TRUE, reload=FALSE)

#' Plot 12: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size
make_figure_A12 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="nlab", epsilon==plot.epsilon, nu==plot.nu, model_name=="BRR",
           n_cal %in% label.values) %>%
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
    plot.file <- sprintf("figures/comparison_delta_marg_const_BRR_nlab_eps%f_nu%f.pdf", plot.epsilon, plot.nu)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(500, 1000, 2000)
label.labels = c("n=500", "n=1000", "n=2000")

make_figure_A12(plot.epsilon=0.1, plot.nu=0.2, save_plots=TRUE, reload=FALSE)

#' Plot 13: Finite sample correction as a function of the calibration set size,
#'         for different number of classes and for different combinations of epsilon and nu
#'

make_figure_A13 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8)) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="Klab", epsilon == plot.epsilon, nu %in% plot.nu.vals, model_name=="BRR") %>%
    filter(K %in% label.values) %>%
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
    plot.file <- sprintf("figures/comparison_delta_marg_const_BRR_Klab_eps%f.pdf", plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(4, 8, 16)
label.labels = c("K=4", "K=8", "K=16")

make_figure_A13(plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8), save_plots=TRUE, reload=FALSE)



#' Plot 14: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size and for different combinations of epsilon and nu
#'         
make_figure_A14 <- function(save_plots=FALSE, reload=FALSE, plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8)) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(plot=="nlab", epsilon == plot.epsilon, nu %in% plot.nu.vals, model_name=="BRR",
           n_cal %in% label.values) %>%
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
    plot.file <- sprintf("figures/comparison_delta_marg_const_BRR_nlab_eps%f.pdf", plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(500, 1000, 2000)
label.labels = c("n=500", "n=1000", "n=2000")

make_figure_A14(plot.epsilon=0.1, plot.nu.vals=c(0.2, 0.8), save_plots=TRUE, reload=FALSE)


# MAIN -------------------------------------------------------------------------
#' Standard --> 1 1
#' Adaptive simplified --> 7 3
#' Adaptive --> 2 0
#' Asymptotic --> 9 5
#' 
#' Adaptive+ simplified --> 10 6
#' Adaptive+ --> 3 2
#' Asymptotic+ --> 4 4
#' 

init_settings <- function(plot.optimistic = FALSE) {
  label.values <<- c("4 classes", "8 classes", "16 classes")
  label.labels <<- c("4 classes", "8 classes", "16 classes")
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized", "Adaptive optimized+")
    method.labels <<- c("Standard", "Adaptive", "Adaptive+")
    color.scale <<- cbPalette[c(1,2,3)]
    shape.scale <<- c(1,0,2)
    linetype.scale <<- c(1,1,1)
    # method.values <<- c("Standard","Adaptive simplified+", "Adaptive optimized+", "Asymptotic+")
    # method.labels <<- c("Standard","Adaptive+ (simplified)", "Adaptive+", "Adaptive+ (asymptotic)")
    # color.scale <<- cbPalette[c(1,8,3,4)]
    # shape.scale <<- c(1,6,2,4)
    # linetype.scale <<- c(1,1,1,1)
  } else {
    method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,7,2,9)]
    shape.scale <<- c(1,3,0,5)
    linetype.scale <<- c(1,1,1,1)
  }
}

load_data <- function(exp.num, from_cluster=TRUE) {
  if(from_cluster) {
    idir <- sprintf("results_hpc/exp%d", exp.num)
  } else {
    idir <- sprintf("results/exp%d", exp.num)
  }        
  ifile.list <- list.files(idir, recursive = FALSE) 
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))    
  summary <- results %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, signal, model_name, contamination, epsilon, nu, estimate, n_train, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

### Experiment 1: UNIFORM -------------------------------------------------------
#' Figure A2
#' Plot marginal coverage (for marginal calibration) as a function of the calibration set size,
#' for different number of classes

make_figure_1 <- function(exp.num=1, plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon=0.1, plot.nu=0.2,
                          plot.optimistic = FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }

  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==20, n_train==10000, signal==1, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, nu==plot.nu) %>%
    filter(n_cal >= 100)
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(K_lab = sprintf("%d classes", K)) %>%
    mutate(K_lab = factor(K_lab, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    #        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
    facet_grid(Key~K_lab, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=8, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}
exp.num <- 1
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.nu <- 0
plot.contamination <- "uniform"

make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

### Experiment 2: Block randomized response model -------------------------------------------------------
#' Figure A3
#' Plot marginal coverage (for marginal calibration) as a function of the calibration set size,
#' for different number of classes

exp.num <- 2
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.nu <- 0
plot.contamination <- "block"

make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


### Experiment 3: 2 level RRB ------------------------------------------------------------
#' Figure A4 and A4.bis
#' Plot marginal coverage (for marginal calibration) as a function of the calibration set
#' size, for different number of classes
exp.num <- 3
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.contamination <- "RRB"
plot.nu <- 0.2

make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

plot.nu <- 0.8
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

### Experiment 4: Uniform, K=4 ------------------------

# init_settings <- function(plot.optimistic = FALSE) {
#   if(plot.optimistic) {
#     method.values <<- c("Standard", "Adaptive optimized+", "Adaptive simplified+", "Asymptotic+")
#     method.labels <<- c("Standard", "Adaptive-o+", "Adaptive-s+", "Asymptotic+")
#   } else {
#     method.values <<- c("Standard", "Adaptive optimized", "Adaptive simplified", "Asymptotic")
#     method.labels <<- c("Standard", "Adaptive-o", "Adaptive-s", "Asymptotic")
#   }
#   label.values <<- c("4 classes", "8 classes", "16 classes")
#   label.labels <<- c("4 classes", "8 classes", "16 classes")
#   cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#   df.dummy <<- tibble(key="Coverage", value=0.95)
#   df.dummy2 <<- tibble(key="Coverage", value=0.5)
#   color.scale <<- cbPalette[c(1,3,7,4)]
#   shape.scale <<- c(1,2,0,3)
#   linetype.scale <<- c(1,1,1,1)
# }

init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  label.values <<- c("4 classes", "8 classes", "16 classes")
  label.labels <<- c("4 classes", "8 classes", "16 classes")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized", "Adaptive optimized+")
    method.labels <<- c("Standard", "Adaptive", "Adaptive+")
    color.scale <<- cbPalette[c(1,2,3)]
    shape.scale <<- c(1,0,2)
    linetype.scale <<- c(1,1,1)
    # method.values <<- c("Standard","Adaptive simplified+", "Adaptive optimized+", "Asymptotic+")
    # method.labels <<- c("Standard","Adaptive+ (simplified)", "Adaptive+", "Adaptive+ (asymptotic)")
    # color.scale <<- cbPalette[c(1,8,3,4)]
    # shape.scale <<- c(1,6,2,4)
    # linetype.scale <<- c(1,1,1,1)
  } else {
    method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,7,2,9)]
    shape.scale <<- c(1,3,0,5)
    linetype.scale <<- c(1,1,1,1)
  }
}

#' Figure 1
#' Plot marginal coverage and size as function of the number of calibration samples.
#' We get two aligned panels.
#' 

make_figure_2 <- function(exp.num=4, plot.alpha=0.1, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon=0.1, plot.nu=0,
                          plot.optimistic = FALSE,
                          slides = FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  if(!slides){
    df <- summary %>%
      filter(data=="synthetic1", num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
             Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
             Method %in% method.values,
             contamination==plot.contamination,
             epsilon==plot.epsilon, nu==plot.nu) %>%
      filter(n_cal >= 500 & n_cal <= 10000)
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
    pp <- df %>%
      mutate(Method = factor(Method, method.values, method.labels)) %>%
      ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.05) +
      facet_wrap(~Key, scales="free_y") +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      scale_x_continuous(trans='log10') +
      xlab("Number of calibration samples") +
      ylab("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    if(save_plots) {
      plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_K%d_eps%f_nu%s_%s_%s_optimistic%s.pdf",
                           exp.num, 10000, plot.K, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=2.25, width=7, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(data == "synthetic1", num_var == 20, n_train == 10000, K==plot.K, signal == 1,
             Guarantee == plot.guarantee, Label == "marginal", model_name == "RFC",
             Alpha == plot.alpha,
             Method %in% method.values,
             contamination == plot.contamination, epsilon == plot.epsilon,
             nu == plot.nu) %>%
      filter(n_cal >= 500 & n_cal <= 10000)

    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
    
    for (i in 1:length(method.values)) {
      current_methods <- method.values[1:i]
      current_labels <- method.labels[1:i]
      
      df_filtered <- df_filt %>%
        filter(Method %in% current_methods) %>%
        mutate(Method = factor(Method, levels = current_methods, labels = current_labels))
      
      {
        df3 = df2 = df_filtered[1:2,]
        df3$n_cal[1] = df2$n_cal[1] = min(df_filtered$n_cal)
        df3$n_cal[2] = df2$n_cal[2] = max(df_filtered$n_cal)
        df2$Mean[1] = 0.8
        df2$Mean[2] = 1.4
        df3$Mean[1] = 1
        df3$Mean[2] = 2
      }
      
      pp <- df_filtered %>%
        ggplot(aes(x = n_cal, y = Mean, color = Method, shape = Method, linetype = Method)) +
        geom_point() +
        geom_line() +
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
        geom_point(data = df2, alpha = 0) +
        geom_point(data = df3, alpha = 0) +
        facet_wrap(~Key, scales = "free_y") +
        geom_hline(data = df.nominal, aes(yintercept = Mean), linetype = "dashed") +
        geom_point(data = df.range, aes(x = n_cal, y = Mean), alpha = 0) +
        scale_color_manual(values = color.scale[1:i]) +
        scale_shape_manual(values = shape.scale[1:i]) +
        scale_linetype_manual(values = linetype.scale[1:i]) +
        scale_x_continuous(trans = 'log10') +
        xlab("Number of calibration samples") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/exp%d_synthetic1_ntrain%d_K%d_eps%f_nu%s_%s_%s_optimistic%s_%d.pdf",
                           exp.num, 10000, plot.K, plot.epsilon, plot.nu, plot.guarantee,
                           plot.contamination, plot.optimistic, i)
      ggsave(file = plot.file, plot = pp, height = 3.5, width = 7, units = "in")
    }
  }
  
}

exp.num = 4
plot.alpha = 0.1
plot.epsilon = 0.1
plot.K = 4

make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.nu=0, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.nu=0, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)
make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.nu=0, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=TRUE)
make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.nu=0, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=TRUE)

#' ---------------------------------------------------------------------------------------------------------------------
#' Figure 2
#' Plot marginal coverage as function of the number of calibration samples, increasing the level of contamination
#' 
 

make_figure_3 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon, plot.nu=0,
                          plot.optimistic=FALSE,
                          slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  if(!slides){
    df <- summary %>%
      filter(data==plot.data, num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
             Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
             Method %in% method.values,
             contamination==plot.contamination,
             nu==plot.nu, epsilon %in% plot.epsilon) %>%
      filter(n_cal >= 500)
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
    pp <- df %>%
      mutate(Method = factor(Method, method.values, method.labels)) %>%
      mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
      #        mutate(Label = factor(Label, label.values, label.labels)) %>%
      ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
      facet_grid(Key~Epsilon, scales="free") +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
      scale_x_continuous(trans='log10') +
      xlab("Number of calibration samples") +
      ylab("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_nu%s_%s_%s_optimistic%s.pdf",
                           exp.num, plot.data, 10000, plot.K, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=3.5, width=7, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(data == plot.data, num_var == 20, n_train == 10000, K==plot.K, signal == 1,
             Guarantee == plot.guarantee, Label == "marginal", model_name == "RFC",
             Alpha == plot.alpha,
             Method %in% method.values,
             contamination == plot.contamination,
             epsilon %in% plot.epsilon, nu == plot.nu) %>%
      filter(n_cal >= 500) %>%
      mutate(Epsilon = sprintf("Contam: %.2f", epsilon))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
    
    for (i in 1:length(method.values)) {
      current_methods <- method.values[1:i]
      current_labels <- method.labels[1:i]
      
      df_filtered <- df_filt %>%
        filter(Method %in% current_methods) %>%
        mutate(Method = factor(Method, levels = current_methods, labels = current_labels)) 
      
      {
        df3 = df2 = df_filtered[1:2,]
        df3$n_cal[1] = df2$n_cal[1] = min(df_filtered$n_cal)
        df3$n_cal[2] = df2$n_cal[2] = max(df_filtered$n_cal)
        df2$Mean[1] = 0.8
        df2$Mean[2] = 1.4
        df3$Mean[1] = 1
        df3$Mean[2] = 3.5
        }
      
      pp <- df_filtered %>%
        ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
        geom_point() +
        geom_line() +
        geom_point(data = df2, alpha = 0) +
        geom_point(data = df3, alpha = 0) +
        facet_grid(Key~Epsilon, scales = "free") +
        geom_hline(data = df.nominal, aes(yintercept = Mean), linetype = "dashed") +
        geom_point(data = df.range, aes(x = n_cal, y = Mean), alpha = 0) +
        scale_color_manual(values = color.scale[1:i]) +
        scale_shape_manual(values = shape.scale[1:i]) +
        scale_linetype_manual(values = linetype.scale[1:i]) +
        scale_x_continuous(trans = 'log10') +
        xlab("Number of calibration samples") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/exp%d_%s_ntrain%d_K%d_nu%s_%s_%s_optimistic%s_%d.pdf",
                           exp.num, plot.data, 10000, plot.K, plot.nu, plot.guarantee,
                           plot.contamination, plot.optimistic, i)
      ggsave(file = plot.file, plot = pp, height = 4, width = 7, units = "in")
      
    }
  }
  
}

exp.num <- 4
plot.alpha <- 0.1
plot.nu <- 0
plot.epsilon <- c(0,0.05,0.1,0.2)
plot.K <- 4
plot.data <- "synthetic1"

make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.epsilon=plot.epsilon, plot.nu=0, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.epsilon=plot.epsilon, plot.nu=0, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.epsilon=plot.epsilon, plot.nu=0, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=TRUE)
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.epsilon=plot.epsilon, plot.nu=0, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=TRUE)

#' Figure A1
#' Plot marginal coverage as function of the strength of label contamination, increasing the number
#' of calibration samples
#' 

make_figure_4 <- function(exp.num, plot.alpha, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.ncal, plot.nu=0,
                          plot.optimistic=FALSE,
                          slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  if(!slides){
    df <- summary %>%
      filter(data=="synthetic1", num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
             Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
             Method %in% method.values,
             contamination==plot.contamination,
             nu==plot.nu, n_cal %in% plot.ncal)
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), epsilon=0.1, Method="Standard")
    pp <- df %>%
      mutate(Method = factor(Method, method.values, method.labels)) %>%
      mutate(Ncal = sprintf("Cal. samples: %d", n_cal)) %>%
      #        mutate(Label = factor(Label, label.values, label.labels)) %>%
      ggplot(aes(x=epsilon, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.01) +
      facet_grid(Key~Ncal, scales="free") +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=epsilon, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
      # scale_x_continuous(trans='log10') +
      xlab("Strength of label contamination") +
      ylab("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/Aexp%d_synthetic1_ntrain%d_K%d_nu%s_%s_%s_optimistic%s.pdf",
                           exp.num, 10000, plot.K, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=3.5, width=7, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(data == "synthetic1", num_var == 20, n_train == 10000, K==plot.K, signal == 1,
             Guarantee == plot.guarantee, Label == "marginal", model_name == "RFC",
             Alpha == plot.alpha,
             Method %in% method.values,
             contamination == plot.contamination,
             n_cal %in% plot.ncal, nu == plot.nu) %>%
      mutate(Ncal = sprintf("Cal. samples: %d", n_cal))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), epsilon=0.1, Method="Standard")
    
    for (i in 1:length(method.values)) {
      current_methods <- method.values[1:i]
      current_labels <- method.labels[1:i]
      
      df_filtered <- df_filt %>%
        filter(Method %in% current_methods) %>%
        mutate(Method = factor(Method, levels = current_methods, labels = current_labels)) 
      
      {
        df3 = df2 = df_filtered[1:2,]
        df3$n_cal[1] = df2$n_cal[1] = min(df_filtered$n_cal)
        df3$n_cal[2] = df2$n_cal[2] = max(df_filtered$n_cal)
        df2$Mean[1] = 0.8
        df2$Mean[2] = 1.4
        df3$Mean[1] = 1
        df3$Mean[2] = 2.4
        }
      
      pp <- df_filtered %>%
        ggplot(aes(x=epsilon, y=Mean, color=Method, shape=Method, linetype=Method)) +
        geom_point() +
        geom_line() +
        geom_point(data = df2, alpha = 0) +
        geom_point(data = df3, alpha = 0) +
        facet_grid(Key~Ncal, scales = "free") +
        geom_hline(data = df.nominal, aes(yintercept = Mean), linetype = "dashed") +
        geom_point(data = df.range, aes(x = epsilon, y = Mean), alpha = 0) +
        scale_color_manual(values = color.scale[1:i]) +
        scale_shape_manual(values = shape.scale[1:i]) +
        scale_linetype_manual(values = linetype.scale[1:i]) +
        # scale_x_continuous(trans = 'log10') +
        xlab("Strength of label contamination") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/Aexp%d_synthetic1_ntrain%d_K%d_nu%s_%s_%s_optimistic%s_%d.pdf",
                           exp.num, 10000, plot.K, plot.nu, plot.guarantee,
                           plot.contamination, plot.optimistic, i)
      ggsave(file = plot.file, plot = pp, height = 4, width = 7, units = "in")
      
    }
  }
  
}

exp.num <- 4
plot.alpha <- 0.1
plot.nu <- 0
plot.ncal <- c(1000,10000,100000)
plot.K <- 4

make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.ncal=plot.ncal, plot.nu=0, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.ncal=plot.ncal, plot.nu=0, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.ncal=plot.ncal, plot.nu=0, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=TRUE)
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.ncal=plot.ncal, plot.nu=0, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
### Experiment 5: RRB, increasing level of contamination, K=4 ------------------------
#' Figure 3
#' Plot marginal coverage as function of the number of calibration samples, increasing the deviation of contamination
#' from a randomized response model
#' 

make_figure_5 <- function(exp.num, plot.alpha, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon, plot.nu=0,
                          plot.optimistic=FALSE,
                          slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  if(!slides){
    df <- summary %>%
      filter(data=="synthetic1", num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
             Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
             Method %in% method.values,
             contamination==plot.contamination,
             epsilon==plot.epsilon, nu %in% plot.nu) %>%
      filter(n_cal >= 500)
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
    
    appender <- function(string) TeX(paste("$\\nu : $", string))  
    
    pp <- df %>%
      mutate(Method = factor(Method, method.values, method.labels)) %>%
      mutate(Nu = sprintf("%.2f", nu)) %>%
      #        mutate(Label = factor(Label, label.values, label.labels)) %>%
      ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
      facet_grid(Key~Nu, scales="free", labeller = labeller(.default=label_parsed,Nu=appender)) +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
      scale_x_continuous(trans='log10') +
      xlab("Number of calibration samples") +
      ylab("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_K%d_eps%s_%s_%s_optimistic%s.pdf",
                           exp.num, 10000, plot.K, plot.epsilon, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=3.5, width=8, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(data == "synthetic1", num_var == 20, n_train == 10000, K==plot.K, signal == 1,
             Guarantee == plot.guarantee, Label == "marginal", model_name == "RFC",
             Alpha == plot.alpha,
             Method %in% method.values,
             contamination == plot.contamination,
             epsilon==plot.epsilon, nu %in% plot.nu) %>%
      filter(n_cal >= 500) %>%
      mutate(Nu = sprintf("%.2f", nu))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
    
    appender <- function(string) TeX(paste("$\\nu : $", string))
    
    for (i in 1:length(method.values)) {
      current_methods <- method.values[1:i]
      current_labels <- method.labels[1:i]
      
      df_filtered <- df_filt %>%
        filter(Method %in% current_methods) %>%
        mutate(Method = factor(Method, levels = current_methods, labels = current_labels)) 
      
      {
        df3 = df2 = df_filtered[1:2,]
        df3$n_cal[1] = df2$n_cal[1] = min(df_filtered$n_cal)
        df3$n_cal[2] = df2$n_cal[2] = max(df_filtered$n_cal)
        df2$Mean[1] = 0.8
        df2$Mean[2] = 1.4
        df3$Mean[1] = 1
        df3$Mean[2] = 2
        }
      
      pp <- df_filtered %>%
        ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
        geom_point() +
        geom_line() +
        geom_point(data = df2, alpha = 0) +
        geom_point(data = df3, alpha = 0) +
        facet_grid(Key~Nu, scales="free", labeller = labeller(.default=label_parsed,Nu=appender)) +
        geom_hline(data = df.nominal, aes(yintercept = Mean), linetype = "dashed") +
        geom_point(data = df.range, aes(x = n_cal, y = Mean), alpha = 0) +
        scale_color_manual(values = color.scale[1:i]) +
        scale_shape_manual(values = shape.scale[1:i]) +
        scale_linetype_manual(values = linetype.scale[1:i]) +
        scale_x_continuous(trans = 'log10') +
        xlab("Number of calibration samples") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/exp%d_synthetic1_ntrain%d_K%d_eps%s_%s_%s_optimistic%s_%d.pdf",
                           exp.num, 10000, plot.K, plot.epsilon, plot.guarantee,
                           plot.contamination, plot.optimistic, i)
      ggsave(file = plot.file, plot = pp, height = 4, width = 7, units = "in")
      
    }
  }
  
}

exp.num <- 5

plot.alpha <- 0.1
plot.nu <- c(0, 0.25, 0.75, 1)
plot.epsilon <- 0.1
plot.K <- 4

make_figure_5(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
make_figure_5(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)
make_figure_5(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=TRUE)
make_figure_5(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=TRUE)

#' ---------------------------------------------------------------------------------------------------------------------
#' Figure 3 bis
#' Plot marginal coverage as function of the deviation from the randomized response model,
#' increasing the number of calibration samples
#' 

make_figure_6 <- function(exp.num, plot.alpha, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="RRB",
                          plot.ncal, plot.epsilon=0.1,
                          plot.optimistic=FALSE,
                          slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  if(!slides){
    df <- summary %>%
      filter(data=="synthetic1", num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
             Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
             Method %in% method.values,
             contamination==plot.contamination,
             epsilon==plot.epsilon, n_cal %in% plot.ncal)
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), nu=0.5, Method="Standard")
    pp <- df %>%
      mutate(Method = factor(Method, method.values, method.labels)) %>%
      mutate(Ncal = sprintf("Cal. samples: %d", n_cal)) %>%
      #        mutate(Label = factor(Label, label.values, label.labels)) %>%
      ggplot(aes(x=nu, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.01) +
      facet_grid(Key~Ncal, scales="free") +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=nu, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
      # scale_x_continuous(trans='log10') +
      xlab(TeX("$\\nu$")) +
      ylab("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/Aexp%d_synthetic1_ntrain%d_K%d_eps%s_%s_%s_optimistic%s.pdf",
                           exp.num, 10000, plot.K, plot.epsilon, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=3.5, width=8, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(data == "synthetic1", num_var == 20, n_train == 10000, K==plot.K, signal == 1,
             Guarantee == plot.guarantee, Label == "marginal", model_name == "RFC",
             Alpha == plot.alpha,
             Method %in% method.values,
             contamination == plot.contamination,
             epsilon == plot.epsilon, n_cal %in% plot.ncal) %>%
      mutate(Ncal = sprintf("Cal. samples: %d", n_cal))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), nu=0.5, Method="Standard")
    
    for (i in 1:length(method.values)) {
      current_methods <- method.values[1:i]
      current_labels <- method.labels[1:i]
      
      df_filtered <- df_filt %>%
        filter(Method %in% current_methods) %>%
        mutate(Method = factor(Method, levels = current_methods, labels = current_labels)) 
      
      {
        df3 = df2 = df_filtered[1:2,]
        df3$n_cal[1] = df2$n_cal[1] = min(df_filtered$n_cal)
        df3$n_cal[2] = df2$n_cal[2] = max(df_filtered$n_cal)
        df2$Mean[1] = 0.8
        df2$Mean[2] = 1.4
        df3$Mean[1] = 1
        df3$Mean[2] = 2.4
        }
      
      pp <- df_filtered %>%
        ggplot(aes(x=nu, y=Mean, color=Method, shape=Method, linetype=Method)) +
        geom_point() +
        geom_line() +
        geom_point(data = df2, alpha = 0) +
        geom_point(data = df3, alpha = 0) +
        facet_grid(Key~Ncal, scales = "free") +
        geom_hline(data = df.nominal, aes(yintercept = Mean), linetype = "dashed") +
        geom_point(data = df.range, aes(x = nu, y = Mean), alpha = 0) +
        scale_color_manual(values = color.scale[1:i]) +
        scale_shape_manual(values = shape.scale[1:i]) +
        scale_linetype_manual(values = linetype.scale[1:i]) +
        # scale_x_continuous(trans = 'log10') +
        xlab(TeX("$\\nu$")) +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/Aexp%d_synthetic1_ntrain%d_K%d_eps%s_%s_%s_optimistic%s_%d.pdf",
                           exp.num, 10000, plot.K, plot.epsilon, plot.guarantee,
                           plot.contamination, plot.optimistic, i)
      ggsave(file = plot.file, plot = pp, height = 4, width = 7, units = "in")
      
    }
  }
  
}

exp.num <- 5
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.ncal <- c(1000,10000,100000)
plot.K <- 4

make_figure_6(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.ncal=plot.ncal, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
make_figure_6(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.ncal=plot.ncal, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)

plot.epsilon <- 0.2
make_figure_6(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.ncal=plot.ncal, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
make_figure_6(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.ncal=plot.ncal, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)


### Experiment 6: Uniform, changing the classifier, K=4 ------------------------
#' Figure Appendix
#' Plot marginal coverage as function of the number of calibration samples, as the classifier changes
#' 


make_figure_7 <- function(exp.num, plot.alpha, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform", plot.model_name,
                          plot.epsilon, plot.nu=0,
                          plot.optimistic=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
           Label=="marginal", model_name %in% plot.model_name, Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, nu==plot.nu) %>%
    filter(n_cal >= 500)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Model = sprintf("Classifier: %s", model_name)) %>%
    #        mutate(Label = factor(Label, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Model, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_K%d_eps%s_%s_%s_optimistic%s_classifiers.pdf",
                         exp.num, 10000, plot.K, plot.epsilon, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=3.5, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
  
}

exp.num <- 6
plot.alpha <- 0.1
plot.nu <- 0
plot.epsilon <- 0.1
plot.K <- 4
plot.model_name <- c("RFC", "SVC", "NN")

make_figure_7(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
              plot.model_name=plot.model_name,
              plot.guarantee="marginal", plot.contamination="uniform",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
make_figure_7(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
              plot.model_name=plot.model_name,
              plot.guarantee="marginal", plot.contamination="uniform",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)



### Experiment 7: Uniform, changing the data distribution, K=4 ------------------------
#' Figure Appendix
#' Plot marginal coverage as function of the number of calibration samples, increasing the level of contamination
#' 

make_figure_8 <- function(exp.num, plot.data, plot.alpha, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon, plot.nu=0,
                          plot.optimistic=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           nu==plot.nu, epsilon %in% plot.epsilon) %>%
    filter(n_cal >= 500)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
    #        mutate(Label = factor(Label, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Epsilon, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num, plot.data, 10000, plot.K, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=3.5, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
    
}

exp.num <- 7
plot.alpha <- 0.1
plot.nu <- 0
plot.epsilon <- c(0,0.05,0.1,0.2)
plot.K <- 4

plot.data <- c("synthetic2")
make_figure_8(exp.num=exp.num, plot.data=plot.data, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.epsilon=plot.epsilon, plot.nu=0, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
make_figure_8(exp.num=exp.num, plot.data=plot.data, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.epsilon=plot.epsilon, plot.nu=0, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

plot.data <- c("synthetic3")
make_figure_8(exp.num=exp.num, plot.data=plot.data, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.epsilon=plot.epsilon, plot.nu=0, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
make_figure_8(exp.num=exp.num, plot.data=plot.data, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="uniform",
              plot.epsilon=plot.epsilon, plot.nu=0, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


### Experiment 8: Uniform, changing the contamination model ------------------------
#' Figure Appendix
#' Plot marginal coverage as function of the strength of label contamination, increasing the number of labels
#' 

make_figure_9 <- function(exp.num, plot.ncal=1000, plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.nu=0.2,
                          plot.optimistic = FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==20, n_train==10000, n_cal==plot.ncal, signal==1,
           Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           nu==plot.nu)
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), epsilon=0.1, Method="Standard")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(K_lab = sprintf("%d classes", K)) %>%
    mutate(K_lab = factor(K_lab, label.values, label.labels)) %>%
    ggplot(aes(x=epsilon, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.01) +
    facet_grid(Key~K_lab, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=epsilon, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    xlab("Strength of label contamination") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_ncal%d_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num, 10000, plot.ncal, plot.nu, plot.guarantee,
                         plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}


exp.num <- 8
plot.alpha <- 0.1
plot.ncal <- 10000
plot.nu <- 0.2

plot.contamination <- "uniform"
make_figure_9(exp.num=exp.num,
              plot.ncal=plot.ncal,
              plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
make_figure_9(exp.num=exp.num,
              plot.ncal=plot.ncal,
              plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

plot.contamination <- "block"
make_figure_9(exp.num=exp.num,
              plot.ncal=plot.ncal,
              plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
make_figure_9(exp.num=exp.num,
              plot.ncal=plot.ncal,
              plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

plot.contamination <- "RRB"
make_figure_9(exp.num=exp.num,
              plot.ncal=plot.ncal,
              plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
make_figure_9(exp.num=exp.num,
              plot.ncal=plot.ncal,
              plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


### Experiment 9: Figures for Conference (and paper) ------------------------------------------

init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized", "Adaptive optimized+")
    method.labels <<- c("Standard", "Adaptive", "Adaptive+")
    color.scale <<- cbPalette[c(1,2,3)]
    shape.scale <<- c(1,0,2)
    linetype.scale <<- c(1,1,1)
    # method.values <<- c("Standard","Adaptive simplified+", "Adaptive optimized+", "Asymptotic+")
    # method.labels <<- c("Standard","Adaptive+ (simplified)", "Adaptive+", "Adaptive+ (asymptotic)")
    # color.scale <<- cbPalette[c(1,8,3,4)]
    # shape.scale <<- c(1,6,2,4)
    # linetype.scale <<- c(1,1,1,1)
  } else {
    # method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized")
    # method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive")
    # color.scale <<- cbPalette[c(1,7,2)]
    # shape.scale <<- c(1,2,0)
    # linetype.scale <<- c(1,1,1)
    method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,7,2,9)]
    shape.scale <<- c(1,3,0,5)
    linetype.scale <<- c(1,1,1,1)
  }
}

make_figure_10 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                              plot.contamination="uniform",
                              plot.epsilon, plot.nu=0,
                              plot.optimistic=FALSE,
                              slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  if(!slides){
    df <- summary %>%
      filter(data==plot.data, num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
             Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
             Method %in% method.values,
             contamination==plot.contamination,
             nu==plot.nu, epsilon %in% plot.epsilon) %>%
      filter(n_cal >= 500)
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
    pp <- df %>%
      mutate(Method = factor(Method, method.values, method.labels)) %>%
      mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
      #        mutate(Label = factor(Label, label.values, label.labels)) %>%
      ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
      facet_grid(Key~Epsilon, scales="free") +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
      scale_x_continuous(trans='log10') +
      xlab("Number of calibration samples") +
      ylab("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_nu%s_%s_%s_optimistic%s.pdf",
                           exp.num, plot.data, 10000, plot.K, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=3.5, width=8, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(data == plot.data, num_var == 20, n_train == 10000, K==plot.K, signal == 1,
             Guarantee == plot.guarantee, Label == "marginal", model_name == "RFC",
             Alpha == plot.alpha,
             Method %in% method.values,
             contamination == plot.contamination,
             epsilon %in% plot.epsilon, nu == plot.nu) %>%
      filter(n_cal >= 500) %>%
      mutate(Epsilon = sprintf("Contam: %.2f", epsilon))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
    
    for (i in 1:length(method.values)) {
      current_methods <- method.values[1:i]
      current_labels <- method.labels[1:i]
      
      df_filtered <- df_filt %>%
        filter(Method %in% current_methods) %>%
        mutate(Method = factor(Method, levels = current_methods, labels = current_labels)) 
      
      if(!plot.optimistic){
        df3 = df2 = df_filtered[1:2,]
        df3$n_cal[1] = df2$n_cal[1] = min(df_filtered$n_cal)
        df3$n_cal[2] = df2$n_cal[2] = max(df_filtered$n_cal)
        df2$Mean[1] = 0.8
        df2$Mean[2] = 1.4
        df3$Mean[1] = 1
        df3$Mean[2] = 2.8
      } else {
        df3 = df2 = df_filtered[1:2,]
        df3$n_cal[1] = df2$n_cal[1] = min(df_filtered$n_cal)
        df3$n_cal[2] = df2$n_cal[2] = max(df_filtered$n_cal)
        df2$Mean[1] = 0.8
        df2$Mean[2] = 1.4
        df3$Mean[1] = 1
        df3$Mean[2] = 2
      } 
      
      pp <- df_filtered %>%
        ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
        geom_point() +
        geom_line() +
        geom_point(data = df2, alpha = 0) +
        geom_point(data = df3, alpha = 0) +
        facet_grid(Key~Epsilon, scales = "free") +
        geom_hline(data = df.nominal, aes(yintercept = Mean), linetype = "dashed") +
        geom_point(data = df.range, aes(x = n_cal, y = Mean), alpha = 0) +
        scale_color_manual(values = color.scale[1:i]) +
        scale_shape_manual(values = shape.scale[1:i]) +
        scale_linetype_manual(values = linetype.scale[1:i]) +
        scale_x_continuous(trans = 'log10') +
        xlab("Number of calibration samples") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/exp%d_%s_ntrain%d_K%d_nu%s_%s_%s_optimistic%s_%d.pdf",
                           exp.num, plot.data, 10000, plot.K, plot.nu, plot.guarantee,
                           plot.contamination, plot.optimistic, i)
      ggsave(file = plot.file, plot = pp, height = 4, width = 6, units = "in")
      
    }
  }
  
}

exp.num <- 9
plot.alpha <- 0.1
plot.nu <- 0.2
plot.epsilon <- c(0,0.05,0.1,0.2)
plot.K <- 4
plot.contamination <- "block"
make_figure_10(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
               plot.contamination=plot.contamination,
               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
make_figure_10(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
               plot.contamination=plot.contamination,
               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)
make_figure_10(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
               plot.contamination=plot.contamination,
               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=TRUE)
make_figure_10(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
               plot.contamination=plot.contamination,
               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=TRUE)

plot.contamination <- "RRB"
make_figure_10(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
                  plot.contamination=plot.contamination,
                  plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
make_figure_10(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
                  plot.contamination=plot.contamination,
                  plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)
make_figure_10(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
                  plot.contamination=plot.contamination,
                  plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=TRUE)
make_figure_10(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
                  plot.contamination=plot.contamination,
                  plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=TRUE)

#'------------------- Plotting the asymptotic method --------------------
init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)")
    color.scale <<- cbPalette[c(1,3,4)]
    shape.scale <<- c(1,2,4)
    linetype.scale <<- c(1,1,1)
  } else {
    method.values <<- c("Standard", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,3,9)]
    shape.scale <<- c(1,0,5)
    linetype.scale <<- c(1,1,1)
  }
}



make_figure_11 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                           plot.contamination="uniform",
                           plot.epsilon, plot.nu=0,
                           plot.optimistic=FALSE,
                           slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  if(!slides){
    df <- summary %>%
      filter(data==plot.data, num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
             Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
             Method %in% method.values,
             contamination==plot.contamination,
             nu==plot.nu, epsilon %in% plot.epsilon) %>%
      filter(n_cal >= 500)
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
    pp <- df %>%
      mutate(Method = factor(Method, method.values, method.labels)) %>%
      mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
      #        mutate(Label = factor(Label, label.values, label.labels)) %>%
      ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
      facet_grid(Key~Epsilon, scales="free") +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
      scale_x_continuous(trans='log10') +
      xlab("Number of calibration samples") +
      ylab("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_nu%s_%s_%s_optimistic%s_asy.pdf",
                           exp.num, plot.data, 10000, plot.K, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=3.5, width=8, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(data == plot.data, num_var == 20, n_train == 10000, K==plot.K, signal == 1,
             Guarantee == plot.guarantee, Label == "marginal", model_name == "RFC",
             Alpha == plot.alpha,
             Method %in% method.values,
             contamination == plot.contamination,
             epsilon %in% plot.epsilon, nu == plot.nu) %>%
      filter(n_cal >= 500) %>%
      mutate(Epsilon = sprintf("Contam: %.2f", epsilon))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
    
    for (i in 1:length(method.values)) {
      current_methods <- method.values[1:i]
      current_labels <- method.labels[1:i]
      
      df_filtered <- df_filt %>%
        filter(Method %in% current_methods) %>%
        mutate(Method = factor(Method, levels = current_methods, labels = current_labels)) 
        
      df3 = df2 = df_filtered[1:2,]
      df3$n_cal[1] = df2$n_cal[1] = min(df_filtered$n_cal)
      df3$n_cal[2] = df2$n_cal[2] = max(df_filtered$n_cal)
      df2$Mean[1] = 0.8
      df2$Mean[2] = 1.4
      df3$Mean[1] = 1
      df3$Mean[2] = 2
      
      pp <- df_filtered %>%
        ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
        geom_point() +
        geom_line() +
        geom_point(data = df2, alpha = 0) +
        geom_point(data = df3, alpha = 0) +
        facet_grid(Key~Epsilon, scales = "free") +
        geom_hline(data = df.nominal, aes(yintercept = Mean), linetype = "dashed") +
        geom_point(data = df.range, aes(x = n_cal, y = Mean), alpha = 0) +
        scale_color_manual(values = color.scale[1:i]) +
        scale_shape_manual(values = shape.scale[1:i]) +
        scale_linetype_manual(values = linetype.scale[1:i]) +
        scale_x_continuous(trans = 'log10') +
        xlab("Number of calibration samples") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/exp%d_%s_ntrain%d_K%d_nu%s_%s_%s_optimistic%s_asy_%d.pdf",
                           exp.num, plot.data, 10000, plot.K, plot.nu, plot.guarantee,
                           plot.contamination, plot.optimistic, i)
      ggsave(file = plot.file, plot = pp, height = 4, width = 6, units = "in")
      
    }
  }
  
}

exp.num <- 9
plot.alpha <- 0.1
plot.nu <- 0.2
plot.epsilon <- c(0,0.05,0.1,0.2)
plot.K <- 4
plot.contamination <- "block"
plot.data <- "synthetic1"
# make_figure_11(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
#                plot.contamination=plot.contamination,
#                plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=FALSE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)
make_figure_11(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
               plot.contamination=plot.contamination,
               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=TRUE)

### Experiment 10: Increasing the number of classes ------------------------------------------------------------
#' Plot marginal coverage (for marginal calibration) as a function of the calibration set
#' size, for different number of classes
#' 
#' 

init_settings <- function(plot.optimistic = FALSE) {
  label.values <<- c("10 classes", "100 classes", "1000 classes")
  label.labels <<- c("10 classes", "100 classes", "1000 classes")
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)")
    color.scale <<- cbPalette[c(1,3,4)]
    shape.scale <<- c(1,2,4)
    linetype.scale <<- c(1,1,1)
    # method.values <<- c("Standard","Adaptive simplified+", "Adaptive optimized+", "Asymptotic+")
    # method.labels <<- c("Standard","Adaptive+ (simplified)", "Adaptive+", "Adaptive+ (asymptotic)")
    # color.scale <<- cbPalette[c(1,8,3,4)]
    # shape.scale <<- c(1,6,2,4)
    # linetype.scale <<- c(1,1,1,1)
  } else {
    method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,7,2,9)]
    shape.scale <<- c(1,3,0,5)
    linetype.scale <<- c(1,1,1,1)
  }
}

make_figure_12 <- function(exp.num=10, plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon=0.1, plot.nu=0.2,
                          plot.optimistic = FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==20, n_train==100000, signal==1, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, nu==plot.nu) %>%
    filter(n_cal >= 100)
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=10000, Method="Standard")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(K_lab = sprintf("%d classes", K)) %>%
    mutate(K_lab = factor(K_lab, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    #        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
    facet_grid(Key~K_lab, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=8, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 10
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.contamination <- "RRB"
plot.nu <- 0.2

make_figure_12(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


### Experiment 101: CIFAR-10 data ------------------------
#' Plot marginal coverage as function of the strength of label contamination, increasing the number of labels
#' 

init_settings <- function(plot.optimistic = FALSE) {
  if(plot.optimistic) {
    # method.values <<- c("Standard", "Adaptive optimized+", "Adaptive simplified+", "Asymptotic+")
    # method.labels <<- c("Standard", "Adaptive-o+", "Adaptive-s+", "Asymptotic+")
    method.values <<- c("Standard", "Adaptive+", "Adaptive optimized+", "Asymptotic+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive-o+", "Asymptotic+")
  } else {
    # method.values <<- c("Standard", "Adaptive optimized", "Adaptive simplified", "Asymptotic")
    # method.labels <<- c("Standard", "Adaptive-o", "Adaptive-s", "Asymptotic")
    method.values <<- c("Standard", "Adaptive", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive", "Adaptive-o", "Asymptotic")
  }
  label.values <<- c("4 classes", "8 classes", "16 classes")
  label.labels <<- c("4 classes", "8 classes", "16 classes")
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  color.scale <<- cbPalette[c(1,3,7,4)]
  shape.scale <<- c(1,2,0,3)
  linetype.scale <<- c(1,1,1,1)
}

# init_settings <- function(plot.optimistic=FALSE) {
#   if(plot.optimistic) {
#     method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+")
#     method.labels <<- c("Standard", "Adaptive-o+", "Asymptotic+")
#   } else {
#     method.values <<- c("Standard", "Adaptive optimized", "Asymptotic")
#     method.labels <<- c("Standard", "Adaptive-o", "Asymptotic")
#   }
#   cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#   df.dummy <- tibble(key="Coverage", value=0.95)
#   df.dummy2 <- tibble(key="Coverage", value=0.5)
#   color.scale <- cbPalette[c(1,3,4)]
#   shape.scale <- c(1,2,3)
#   linetype.scale <- c(1,1,1)
# }

load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
  ifile.list <- list.files(idir)
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  
  summary <- results %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, K, n_cal, n_test, epsilon_n_clean, epsilon_n_corr, estimate, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
  
  return(summary)
}

make_figure_101 <- function(exp.num, plot.alpha=0.1, plot.K, plot.estimate="rho-epsilon-point",
                           plot.guarantee="marginal",
                           plot.optimistic=FALSE, save_plots=FALSE, reload=TRUE,
                           slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic=plot.optimistic)
  
  if(!slides){
    
    df <- summary %>%
      filter(Alpha==plot.alpha, K==plot.K, estimate==plot.estimate, Guarantee==plot.guarantee, Label=="marginal",
             Method %in% method.values, n_cal %in% c(500,1500,4500,9500))  %>%
      mutate(Method = factor(Method, method.values, method.labels))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.88,0.94), n_cal=1000, Method="Standard")
    
    pp <- df %>%
      ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.1) +
      facet_wrap(.~Key, scales="free") +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      scale_x_continuous(trans='log10', limits=c(500,10000)) +
      xlab("Number of calibration samples") +
      ylab("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    if(save_plots) {
      plot.file <- sprintf("figures/cifar10_%s_optimistic%s_%s.pdf", plot.guarantee, plot.optimistic, plot.estimate)
      ggsave(file=plot.file, height=2.25, width=6.5, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(Alpha==plot.alpha, K==plot.K, estimate==plot.estimate, Guarantee==plot.guarantee, Label=="marginal",
             Method %in% method.values, n_cal %in% c(500,1500,4500,9500))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.88,0.94), n_cal=5000, Method="Standard")
    
    for (i in 1:length(method.values)) {
      current_methods <- method.values[1:i]
      current_labels <- method.labels[1:i]
      
      df_filtered <- df_filt %>%
        filter(Method %in% current_methods) %>%
        mutate(Method = factor(Method, levels = current_methods, labels = current_labels))
      
      {
        df3 = df2 = df_filtered[1:2,]
        df3$n_cal[1] = df2$n_cal[1] = min(df_filtered$n_cal)
        df3$n_cal[2] = df2$n_cal[2] = max(df_filtered$n_cal)
        df2$Mean[1] = 0.88
        df2$Mean[2] = 1.05
        df3$Mean[1] = 0.94
        df3$Mean[2] = 1.25
        }
      
      pp <- df_filtered %>%
        ggplot(aes(x = n_cal, y = Mean, color = Method, shape = Method, linetype = Method)) +
        geom_point() +
        geom_line() +
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
        geom_point(data = df2, alpha = 0) +
        geom_point(data = df3, alpha = 0) +
        facet_wrap(~Key, scales = "free_y") +
        geom_hline(data = df.nominal, aes(yintercept = Mean), linetype = "dashed") +
        geom_point(data = df.range, aes(x = n_cal, y = Mean), alpha = 0) +
        scale_color_manual(values = color.scale[1:i]) +
        scale_shape_manual(values = shape.scale[1:i]) +
        scale_linetype_manual(values = linetype.scale[1:i]) +
        scale_x_continuous(trans='log10', limits=c(500,10000)) +
        xlab("Number of calibration samples") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/cifar10_%s_optimistic%s_%s_%d.pdf",
                           plot.guarantee, plot.optimistic, plot.estimate, i)
      ggsave(file = plot.file, plot = pp, height = 3.5, width = 7, units = "in")
    }
  }
  
  
}


exp.num <- 101
plot.alpha <- 0.1
plot.K <- 10

make_figure_101(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.estimate="rho-epsilon-point", plot.guarantee="marginal",
               plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)

make_figure_101(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.estimate="rho", plot.guarantee="marginal",
               plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)


## make_figure_101(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal",
##               plot.optimistic=FALSE, save_plots=TRUE, reload=TRUE)
## make_figure_101(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal",
##                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)
## make_figure_101(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal",
##                plot.optimistic=FALSE, save_plots=TRUE, reload=TRUE, slides=TRUE)
## make_figure_101(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal",
##                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE, slides=TRUE)


### Experiment 201: BigEarthNet data ------------------------
#' Plot marginal coverage as function of the strength of label contamination, increasing the number of labels
#' 

init_settings <- function(plot.optimistic = FALSE) {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  if(plot.optimistic) {
    # method.values <<- c("Standard", "Adaptive optimized+", "Adaptive simplified+", "Asymptotic+")
    # method.labels <<- c("Standard", "Adaptive-o+", "Adaptive-s+", "Asymptotic+")
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)")
    color.scale <<- cbPalette[c(1,3,4)]
    shape.scale <<- c(1,2,4)
    linetype.scale <<- c(1,1,1)
  } else {
    # method.values <<- c("Standard", "Adaptive optimized", "Adaptive simplified", "Asymptotic")
    # method.labels <<- c("Standard", "Adaptive-o", "Adaptive-s", "Asymptotic")
    method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,7,2,9)]
    shape.scale <<- c(1,3,0,5)
    linetype.scale <<- c(1,1,1,1)
  }
  label.values <<- c("4 classes", "8 classes", "16 classes")
  label.labels <<- c("4 classes", "8 classes", "16 classes")
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
}


load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
  ifile.list <- list.files(idir)
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  
  summary <- results %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, K, n_cal, n_test, epsilon_n_clean, epsilon_n_corr,
             estimate, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
  
  return(summary)
}

make_figure_201 <- function(exp.num, plot.alpha=0.1, plot.K, plot.estimate="rho",
                            plot.guarantee="marginal",
                            plot.optimistic=FALSE, save_plots=FALSE, reload=TRUE,
                            slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic=plot.optimistic)
  
  if(!slides){
    
    df <- summary %>%
      filter(Alpha==plot.alpha, K==plot.K, estimate==plot.estimate, Guarantee==plot.guarantee, Label=="marginal",
             Method %in% method.values, n_cal %in% c(500,1500,2500))  %>%
      mutate(Method = factor(Method, method.values, method.labels))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.88,0.94), n_cal=1500, Method="Standard")
    
    pp <- df %>%
      ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.1) +
      facet_wrap(.~Key, scales="free") +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      scale_x_continuous(trans='log10', limits=c(500,3000)) +
      xlab("Number of calibration samples") +
      ylab("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    if(save_plots) {
      plot.file <- sprintf("figures/bigearthnet_K%d_%s_optimistic%s_%s.pdf", plot.K, plot.guarantee, plot.optimistic, plot.estimate)
      ggsave(file=plot.file, height=2.25, width=6.5, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(Alpha==plot.alpha, K==plot.K, estimate==plot.estimate, Guarantee==plot.guarantee, Label=="marginal",
             Method %in% method.values, n_cal %in% c(500,1500,2500))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.88,0.94), n_cal=1500, Method="Standard")
    
    for (i in 1:length(method.values)) {
      current_methods <- method.values[1:i]
      current_labels <- method.labels[1:i]
      
      df_filtered <- df_filt %>%
        filter(Method %in% current_methods) %>%
        mutate(Method = factor(Method, levels = current_methods, labels = current_labels))
      
      {
        df3 = df2 = df_filtered[1:2,]
        df3$n_cal[1] = df2$n_cal[1] = min(df_filtered$n_cal)
        df3$n_cal[2] = df2$n_cal[2] = max(df_filtered$n_cal)
        df2$Mean[1] = 0.88
        df2$Mean[2] = 1.05
        df3$Mean[1] = 0.94
        df3$Mean[2] = 1.25
        }
      
      pp <- df_filtered %>%
        ggplot(aes(x = n_cal, y = Mean, color = Method, shape = Method, linetype = Method)) +
        geom_point() +
        geom_line() +
        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
        geom_point(data = df2, alpha = 0) +
        geom_point(data = df3, alpha = 0) +
        facet_wrap(~Key, scales = "free_y") +
        geom_hline(data = df.nominal, aes(yintercept = Mean), linetype = "dashed") +
        geom_point(data = df.range, aes(x = n_cal, y = Mean), alpha = 0) +
        scale_color_manual(values = color.scale[1:i]) +
        scale_shape_manual(values = shape.scale[1:i]) +
        scale_linetype_manual(values = linetype.scale[1:i]) +
        scale_x_continuous(trans='log10', limits=c(500,3000)) +
        xlab("Number of calibration samples") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/bigearthnet_K%d_%s_optimistic%s_%s_%d.pdf",
                           plot.K, plot.guarantee, plot.optimistic, plot.estimate, i)
      ggsave(file = plot.file, plot = pp, height = 3.5, width = 7, units = "in")
    }
  }
  
  
}


exp.num <- 201
plot.alpha <- 0.1
plot.K <- 6
plot.epsilon <- 0.017

make_figure_201(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.epsilon=plot.epsilon,
                plot.estimate="rho", plot.guarantee="marginal",
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)


