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
### Experiment 1: UNIFORM -------------------------------------------------------
load_data <- function(exp.num) {
  idir <- sprintf("results/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE)
  
  # Filtra i file per escludere quelli nella sottocartella simplified_methods
  ifile.list <- ifile.list[!grepl("simplified_methods", ifile.list)]
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  
  
  summary <- results %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, signal, model_name, contamination, epsilon, nu, estimate, n_train, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
  
  return(summary)
}

exp.num <- 1
summary <- load_data(exp.num)

# method.values = c("Standard", "Adaptive", "Adaptive optimized", "Adaptive simplified")
# method.labels = c("Standard", "Adaptive(Old)", "Adaptive-opt", "Adaptive-simpl")
method.values = c("Standard", "Adaptive optimized", "Adaptive simplified")
method.labels = c("Standard", "Adaptive-opt", "Adaptive-simpl")
label.values = c("4 classes", "8 classes", "16 classes")
label.labels = c("4 classes", "8 classes", "16 classes")
cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.5)
# color.scale <- cbPalette[c(1,8,3,7)]
# shape.scale <- c(1,3,2,0)
# linetype.scale <- c(1,1,1,1)
color.scale <- cbPalette[c(1,3,7)]
shape.scale <- c(1,2,0)
linetype.scale <- c(1,1,1)

#' --------------------------------------------------
#' Plot marginal coverage (for marginal calibration) as a function of the calibration set size,
#' for different number of classes

make_figure_1 <- function(plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon=0.1, plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==10, n_train==1000, signal==1, Guarantee==plot.guarantee,
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
    plot.file <- sprintf("figures/synthetic1_ntrain%d_eps%f_nu%s_%s_%s.pdf", 10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.nu <- "none"
plot.contamination <- "uniform"

make_figure_1(plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, reload=FALSE)

### Experiment 1: Block randomized response model -------------------------------------------------------
load_data <- function(exp.num) {
  idir <- sprintf("results/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE)
  
  # Filtra i file per escludere quelli nella sottocartella simplified_methods
  ifile.list <- ifile.list[!grepl("simplified_methods", ifile.list)]
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  
  
  summary <- results %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, signal, model_name, contamination, epsilon, nu, estimate, n_train, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
  
  return(summary)
}

exp.num <- 1
summary <- load_data(exp.num)

# method.values = c("Standard", "Adaptive", "Adaptive optimized", "Adaptive simplified")
# method.labels = c("Standard", "Adaptive(Old)", "Adaptive-opt", "Adaptive-simpl")
method.values = c("Standard", "Adaptive optimized", "Adaptive simplified")
method.labels = c("Standard", "Adaptive-opt", "Adaptive-simpl")
label.values = c("4 classes", "8 classes", "16 classes")
label.labels = c("4 classes", "8 classes", "16 classes")
cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.5)
# color.scale <- cbPalette[c(1,8,3,7)]
# shape.scale <- c(1,3,2,0)
# linetype.scale <- c(1,1,1,1)
color.scale <- cbPalette[c(1,3,7)]
shape.scale <- c(1,2,0)
linetype.scale <- c(1,1,1)

#' --------------------------------------------------
#' Plot marginal coverage (for marginal calibration) as a function of the calibration set size,
#' for different number of classes

make_figure_2 <- function(plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon=0.1, plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==10, n_train==1000, signal==1, Guarantee==plot.guarantee,
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
    plot.file <- sprintf("figures/synthetic1_ntrain%d_eps%f_nu%s_%s_%s.pdf", 10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.nu <- "none"
plot.contamination <- "block"

make_figure_2(plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, reload=FALSE)


### Experiment 1: 2 level RRB ------------------------------------------------------------

load_data <- function(exp.num) {
  idir <- sprintf("results/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE)
  
  # Filtra i file per escludere quelli nella sottocartella simplified_methods
  ifile.list <- ifile.list[!grepl("simplified_methods", ifile.list)]
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  
  
  summary <- results %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, signal, model_name, contamination, epsilon, nu, estimate, n_train, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
  
  return(summary)
}

exp.num <- 1
summary <- load_data(exp.num)

# method.values = c("Standard", "Adaptive", "Adaptive optimized", "Adaptive simplified")
# method.labels = c("Standard", "Adaptive(Old)", "Adaptive-opt", "Adaptive-simpl")
method.values = c("Standard", "Adaptive optimized", "Adaptive simplified")
method.labels = c("Standard", "Adaptive-opt", "Adaptive-simpl")
label.values = c("4 classes", "8 classes", "16 classes")
label.labels = c("4 classes", "8 classes", "16 classes")
cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df.dummy <- tibble(key="Coverage", value=0.95)
df.dummy2 <- tibble(key="Coverage", value=0.5)
# color.scale <- cbPalette[c(1,8,3,7)]
# shape.scale <- c(1,3,2,0)
# linetype.scale <- c(1,1,1,1)
color.scale <- cbPalette[c(1,3,7)]
shape.scale <- c(1,2,0)
linetype.scale <- c(1,1,1)

#' --------------------------------------------------
#' Plot marginal coverage (for marginal calibration) as a function of the calibration set
#' size, for different number of classes

make_figure_3 <- function(plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon=0.1, plot.nu=0.2) {
  if(reload) {
    summary <- load_data(1)
  }
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==10, n_train==1000, signal==1, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           contamination==plot.contamination, Method %in% method.values,
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
    plot.file <- sprintf("figures/synthetic1_ntrain%d_eps%f_nu%f_%s_%s.pdf", 10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.nu <- 0.2
plot.contamination <- "RRB"

make_figure_3(plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, reload=FALSE)



