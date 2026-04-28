options(width = 300)

library(tidyverse)
library(latex2exp)
library(RColorBrewer)


### Experiment 1: Impact of Label contamination strength -----------------------
#' Figure 1, Figure A9 and Figure A10
#' Plot marginal coverage as function of the number of calibration samples, increasing the strength
#' of the label contamination
#' 
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
  } else {
    method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,7,2,9)]
    shape.scale <<- c(1,3,0,5)
    linetype.scale <<- c(1,1,1,1)
  }
}

make_figure_1 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
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
             nu==plot.nu, epsilon %in% plot.epsilon)
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.875,1), n_cal=1000, Method="Standard")
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
      theme(text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            plot.margin = margin(5, 1, 1, -10))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_nu%s_%s_%s_optimistic%s.pdf",
                           exp.num, plot.data, 10000, plot.K, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=3.2, width=9, units="in")
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

exp.num <- 1
plot.alpha <- 0.1
plot.nu <- 0.2
plot.epsilon <- c(0,0.05,0.1,0.2)
plot.K <- 4
plot.data <- "synthetic1"

plot.contamination <- "RRB"
## Figure 1
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
              plot.guarantee="marginal",
              plot.contamination=plot.contamination, plot.epsilon=plot.epsilon, plot.nu=plot.nu,
              save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)

# Optimistic counterpart (not shown in paper)
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)

plot.contamination <- "uniform"
## Figure A9
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)


plot.contamination <- "block"
## Figure A10
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)

make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)

#' ---------------------------------------------------------------------------------------------------------------------
### Experiment 2: Impact of Label contamination model ------------------------
#' Figure 2 and Figure A11
#' Plot marginal coverage as function of the number of calibration samples, increasing the deviation of contamination
#' from a randomized response model
#' 

make_figure_2 <- function(exp.num, plot.alpha, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
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
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.875,0.975), n_cal=1000, Method="Standard")
    
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
      theme(text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            plot.margin = margin(5, 5, 1, -10))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_K%d_eps%s_%s_%s_optimistic%s.pdf",
                           exp.num, 10000, plot.K, plot.epsilon, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=3.2, width=9, units="in")
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
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.85,1), n_cal=1000, Method="Standard")
    
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

exp.num <- 2

plot.alpha <- 0.1
plot.nu <- c(0, 0.25, 0.75, 1)
plot.epsilon <- 0.1
plot.K <- 4

## Figure 2
make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
## Optimistic counterpart (not shown in paper)
# make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)


# # For slides
# make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=TRUE)
# # Optimistic counterpart (not shown in paper)
# make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=TRUE)


make_figure_2A <- function(exp.num, plot.alpha, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
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
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.85,1), nu=0.5, Method="Standard")
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
      theme(text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            plot.margin = margin(5, 5, 1, -10))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/Aexp%d_synthetic1_ntrain%d_K%d_eps%s_%s_%s_optimistic%s.pdf",
                           exp.num, 10000, plot.K, plot.epsilon, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=4, width=7, units="in")
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

exp.num <- 2
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.ncal <- c(1000,10000,100000)
plot.K <- 4

## Figure A11
make_figure_2A(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.ncal=plot.ncal, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
# Optimistic counterpart (not shown in paper)
# make_figure_2A(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
#               plot.ncal=plot.ncal, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)

#' ---------------------------------------------------------------------------------------------------------------------
### Experiment 3: Impact of number of classes, RRB -----------------------------
#' Figure A12 and Figure A13
#' Plot marginal coverage as function of the number of calibration samples, increasing the number of classes
#' 

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
  } else {
    method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,7,2,9)]
    shape.scale <<- c(1,3,0,5)
    linetype.scale <<- c(1,1,1,1)
  }
}

make_figure_3 <- function(exp.num=1, plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
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
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.85,1), n_cal=1000, Method="Standard")
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
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 5, 1, -10))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}


exp.num <- 3
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.contamination <- "RRB"

## Figure A12
plot.nu <- 0.2
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# Optimistic counterpart (not shown in paper)
# make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

## Figure A13
plot.nu <- 0.8
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
## Optimistic counterpart (not shown in paper)
# make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
### Experiment 4: Impact of number of classes, Uniform and Block ------------------------
#' Figure A14 and Figure A15
#' Plot marginal coverage as function of the number of calibration samples, increasing the number of classes
#' 

exp.num <- 4
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.nu <- 0

## Figure A14: Uniform
plot.contamination <- "uniform"
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# Optimistic counterpart (not shown in paper)
# make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

## Figure A15: Block
plot.contamination <- "block"
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# Optimistic counterpart (not shown in paper)
# make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)



### Experiment 6: Impact of higher K ------------------------------------------------

init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  label.values <<- c("10 classes", "20 classes", "50 classes")
  label.labels <<- c("10 classes", "20 classes", "50 classes")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized", "Adaptive optimized+")
    method.labels <<- c("Standard", "Adaptive", "Adaptive+")
    color.scale <<- cbPalette[c(1,2,3)]
    shape.scale <<- c(1,0,2)
    linetype.scale <<- c(1,1,1)
  } else {
    method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,7,2,9)]
    shape.scale <<- c(1,3,0,5)
    linetype.scale <<- c(1,1,1,1)
  }
}

make_figure_4 <- function(exp.num=1, plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
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
           epsilon==plot.epsilon, nu==plot.nu,
           K %in% c(10,20,50)) %>%
    filter(n_cal >= 500)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.875,0.975), n_cal=1000, Method="Standard")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(K_lab = sprintf("%d classes", K)) %>%
    mutate(K_lab = factor(K_lab, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    #        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
    facet_wrap(Key~K_lab, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    #ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))

  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=3, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 6
plot.alpha <- 0.1
plot.nu <- 0.8

plot.contamination <- "uniform"
plot.epsilon <- 0.1
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=FALSE, plot.optimistic=FALSE, reload=TRUE)

plot.epsilon <- 0.2
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)


plot.contamination <- "block"
plot.epsilon <- 0.1
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# # Optimistic counterpart (not shown in paper)
# make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

plot.epsilon <- 0.2
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# # Optimistic counterpart (not shown in paper)
# make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)
# 

plot.contamination <- "RRB"
plot.epsilon <- 0.1
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# # Optimistic counterpart (not shown in paper)
# make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

plot.epsilon <- 0.2
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)

# Figure 3 (extended visualization)
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


#### Horizontal visualization

init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  label.values <<- c("10 classes", "50 classes")
  label.labels <<- c("10 classes", "50 classes")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized", "Adaptive optimized+")
    method.labels <<- c("Standard", "Adaptive", "Adaptive+")
    color.scale <<- cbPalette[c(1,2,3)]
    shape.scale <<- c(1,0,2)
    linetype.scale <<- c(1,1,1)
  } else {
    method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,7,2,9)]
    shape.scale <<- c(1,3,0,5)
    linetype.scale <<- c(1,1,1,1)
  }
}

make_figure_4_horizontal <- function(exp.num, plot.alpha=0.1, plot.guarantee="marginal",
                                     plot.contamination=plot.contamination,
                                     plot.epsilon=plot.epsilon, plot.nu=plot.nu,
                                     plot.optimistic=FALSE, save_plots=FALSE, reload=TRUE) {
  
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic=plot.optimistic)
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==20, n_train==10000, signal==1, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, nu==plot.nu,
           K %in% c(10,50)) %>%
    filter(n_cal >= 500)
  
  df_10 <- df %>% filter(K == 10)
  df_50 <- df %>% filter(K == 50)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.875,0.975), n_cal=1000, Method="Standard")
  
  # 10 Classes panel
  pp_10 <- df_10 %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(K_lab = sprintf("%d classes", K)) %>%
    mutate(K_lab = factor(K_lab, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    #        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
    facet_wrap(K_lab~Key, scales="free") +
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
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  # 50 classes Panel
  pp_50 <- df_50 %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(K_lab = sprintf("%d classes", K)) %>%
    mutate(K_lab = factor(K_lab, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    #        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
    facet_wrap(K_lab~Key, scales="free") +
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
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  # Remove individual x-axis labels from the plots
  pp_10 <- pp_10 + theme(legend.position = "none")
  pp_50 <- pp_50 + theme(legend.position = "none")
  
  # Combine the plots, add a shared legend, and a single x-axis label
  pp <- pp_10 + pp_50 + 
    plot_layout(ncol = 2, guides = "collect") & # Combine guides (legend)
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.margin = margin(2, 2, -5, -10)
    )
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s_horizontal.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=3, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
  
}

exp.num <- 6
plot.alpha <- 0.1
plot.nu <- 0.8

## Figure 3 (horizontal visualization)
plot.contamination <- "RRB"
plot.epsilon <- 0.2
make_figure_4_horizontal(exp.num=exp.num, plot.alpha=plot.alpha,
                         plot.guarantee="marginal", plot.contamination=plot.contamination,
                         plot.epsilon=plot.epsilon, plot.nu=plot.nu,
                         save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


### Experiment 5: The advantage of a method for marginal coverage ------------------------------------------------
#' Contamination model: RRB with nu = 0.2 -------------------------------------
#' Plot marginal coverage as function of the number of calibration samples, increasing the number of classes
#' 

init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  label.values <<- c("10 classes", "20 classes", "50 classes")
  label.labels <<- c("10 classes", "20 classes", "50 classes")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+", "Label conditional+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)", "Adaptive+ (label-cond)")
    color.scale <<- cbPalette[c(1,3,4,8)]
    shape.scale <<- c(1,2,4,7)
    linetype.scale <<- c(1,1,1,1)
  }
}

make_figure_5 <- function(exp.num=1, plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
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
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.85,1), n_cal=1000, Method="Standard")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(K_lab = sprintf("%d classes", K)) %>%
    mutate(K_lab = factor(K_lab, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    #        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
    facet_wrap(Key~K_lab, scales="free_y") +
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
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 5, 1, -10))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=5, width=7.5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 5
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.contamination <- "RRB"
plot.nu <- 0.2

# Not shown in paper, shows similar results to those of experiment 7 (Figure 4)
make_figure_5(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


### Experiment 7: The advantage of a method for marginal coverage ------------------------------------------------
#' Contamination model: RRB with nu = 0.2 -------------------------------------
#' Plot marginal coverage as function of the number of calibration samples, increasing the number of classes
#' The separability of the classes increases as K increases, so that the average size of the
#' prediction sets remains stable across experiments with different K

load_data <- function(exp.num, from_cluster=TRUE) {
  if(from_cluster) {
    idir <- sprintf("results_hpc/exp%d", exp.num)
  } else {
    idir <- sprintf("results/exp%d", exp.num, plot.signal)
  }        
  ifile.list <- list.files(idir, recursive = FALSE) 
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))    
  summary <- results %>%
    filter(signal==plot.signal) %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, model_name, contamination, epsilon, nu, estimate, n_train, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  label.values <<- c("10 classes", "20 classes", "50 classes")
  label.labels <<- c("10 classes", "20 classes", "50 classes")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+", "Label conditional+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)", "Adaptive+ (label-cond)")
    color.scale <<- cbPalette[c(1,3,4,8)]
    shape.scale <<- c(1,2,4,7)
    linetype.scale <<- c(1,1,1,1)
  }
}

make_figure_7 <- function(exp.num=1, plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon=0.1, plot.nu=0.2, plot.signal=5,
                          plot.optimistic = FALSE) {
  if(reload) {
    summary <- load_data(exp.num, plot.signal)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==20, n_train==10000, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, nu==plot.nu, n_cal>=500)
    
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.85,1), n_cal=1000, Method="Standard")
  
  pp <- df %>%
    # mutate(Method = factor(Method, method.values, method.labels),
    #        Mean = ifelse(Key == "Size", log10(Mean), Mean),
    #        SE = ifelse(Key == "Size", log10(Mean + SE) - log10(Mean), SE)) %>%
    mutate(Method = factor(Method, method.values, method.labels),
           Mean = ifelse( (K==50) & (Mean>8), NA, Mean),
           Mean = ifelse( (K==20) & (Mean>2), NA, Mean)) %>%
    mutate(K_lab = factor(sprintf("%d classes", K), 
                          levels = sprintf("%d classes", c(10, 20, 50)), 
                          labels = c("10 classes", "20 classes", "50 classes"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    #        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
    #facet_grid(Key~K_lab, scales="free_y", labeller=custom_labeller)+
    facet_wrap(Key~K_lab, scales="free_y")+
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10')+
    xlab("Number of calibration samples") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 5, 1, -10))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=5, width=7.5, units="in")
    #ggsave(file=plot.file, height=4, width=7.5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}


exp.num <- 7
plot.alpha <- 0.1
plot.epsilon <- 0.05
plot.contamination <- "RRB"
plot.nu <- 0.2
plot.signal <- 5

## Figure 4 (extended visualization)
make_figure_7(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.signal=plot.signal, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


#### Horizontal visualization

init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  label.values <<- c("10 classes", "50 classes")
  label.labels <<- c("10 classes", "50 classes")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+", "Label conditional+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)", "Adaptive+ (label-cond)")
    color.scale <<- cbPalette[c(1,3,4,8)]
    shape.scale <<- c(1,2,4,7)
    linetype.scale <<- c(1,1,1,1)
  }
}

make_figure_7_horizontal <- function(exp.num, plot.alpha=0.1, plot.guarantee="marginal",
                                     plot.contamination=plot.contamination,
                                     plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.signal=5,
                                     plot.optimistic=FALSE, save_plots=FALSE, reload=TRUE) {
  
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic=plot.optimistic)
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==20, n_train==10000, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, nu==plot.nu, K %in% c(10,50), n_cal>=500)
  
  df_10 <- df %>% filter(K == 10)
  df_50 <- df %>% filter(K == 50)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.875,0.975), n_cal=1000, Method="Standard")
  
  # 10 Classes panel
  pp_10 <- df_10 %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(K_lab = sprintf("%d classes", K)) %>%
    mutate(K_lab = factor(K_lab, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    facet_wrap(K_lab~Key, scales="free_y") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    scale_y_continuous(trans='log10')+
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  # 50 classes Panel
  pp_50 <- df_50 %>%
    mutate(Method = factor(Method, method.values, method.labels),
           Mean = ifelse(Mean>8, NA, Mean)) %>%
    mutate(K_lab = sprintf("%d classes", K)) %>%
    mutate(K_lab = factor(K_lab, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    facet_wrap(K_lab~Key, scales="free_y") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    scale_y_continuous(trans='log10')+
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  # Remove individual x-axis labels from the plots
  pp_10 <- pp_10 + theme(legend.position = "none")
  pp_50 <- pp_50 + theme(legend.position = "none")
  
  # Combine the plots, add a shared legend, and a single x-axis label
  pp <- pp_10 + pp_50 + 
    plot_layout(ncol = 2, guides = "collect") & # Combine guides (legend)
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.margin = margin(2, 2, -5, -10)
    )
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s_horizontal.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=3, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
  
}


exp.num <- 7
plot.alpha <- 0.1
plot.epsilon <- 0.05
plot.contamination <- "RRB"
plot.nu <- 0.2
plot.signal <- 5

## Figure 4 (horizontal visualization)
make_figure_7_horizontal(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
                         plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.signal=plot.signal,
                         save_plots=TRUE, plot.optimistic=TRUE, reload=FALSE)

### Experiment 301: RBB, Increase in the class-imbalance -------------------

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
    group_by(data, num_var, K, signal, model_name, contamination, epsilon, nu, imb, estimate, n_train, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
    #summarise(Mean=mean(Value, na.rm=TRUE), N=sum(!is.na(Value)), SE=2*sd(Value, na.rm=TRUE)/sqrt(N))
  return(summary)
}


init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+", "Label conditional+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)", "Adaptive+ (label-cond)")
    color.scale <<- cbPalette[c(1,3,4,8)]
    shape.scale <<- c(1,2,4,7)
    linetype.scale <<- c(1,1,1,1)
  }
}

make_figure_301 <- function(exp.num, plot.alpha, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon, plot.nu=0,
                          imb.values,
                          plot.data="synthetic4",
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
           epsilon==plot.epsilon, nu %in% plot.nu, imb %in% imb.values) %>%
    filter(n_cal >= 500)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.85,1), n_cal=1000, Method="Standard")
  
  appender <- function(string) TeX(paste("$\\mu : $", string))  
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Nu = sprintf("%.2f", imb)) %>%
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
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 5, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_eps%s_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num, plot.data, 10000, plot.K, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=7.5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
  
}


exp.num <- 301
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.contamination <- "RRB"
plot.nu <- 0.2
plot.K <- 4
plot.data <- "synthetic4"
imb.values <- c(0,0.5,1)

## Figure A16
make_figure_301(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, imb.values=imb.values, plot.data=plot.data, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

## Plot of the label-wise performances (not shown in paper)
#'
#'

make_figure_302 <- function(exp.num=exp.num, plot.data="synthetic5", plot.alpha=0.1,
                            plot.K=4, plot.epsilon=0.1, plot.nu=0.2,
                            plot.guarantee="marginal", plot.contamination="uniform",
                            plot.estimate="none", plot.imb=1,
                            plot.optimistic=TRUE, save_plots=FALSE, reload=TRUE){
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic=plot.optimistic)
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train==10000, K==plot.K, estimate==plot.estimate, Guarantee==plot.guarantee,
           model_name=="RFC", Alpha==plot.alpha, epsilon==plot.epsilon, nu==plot.nu, contamination==plot.contamination,
           Method %in% method.values, Label %in% label.values, imb==plot.imb)
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  #df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.88,0.95), n_cal=1000, Method="Standard")
  df.range2 <- tibble(Key=c("Size","Size"), Mean=c(1,2), n_cal=1000, Method="Standard")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Label = factor(Label, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    #        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
    facet_grid(Key~Label, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    geom_point(data=df.range2, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_eps%s_nu%s_mu%s_%s_%s_optimistic%s_labelwise.pdf",
                         exp.num, plot.data, 10000, plot.K, plot.epsilon, plot.nu, plot.imb, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}


exp.num <- 301
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.contamination <- "RRB"
plot.nu <- 0.2
plot.K <- 4
plot.data <- "synthetic4"
plot.imb=2

label.values <<- 0:3
label.labels <<- paste("Label", 1:4, sep=" ")

make_figure_302(exp.num=exp.num, plot.data=plot.data, plot.alpha=plot.alpha,
                plot.K=plot.K, plot.epsilon=plot.epsilon, plot.nu=plot.nu,
                plot.guarantee="marginal", plot.contamination=plot.contamination,
                plot.estimate="none", plot.imb=plot.imb,
                plot.optimistic=TRUE, save_plots=FALSE, reload=TRUE)


### Experiment 101: CIFAR-10H data ------------------------
#' Figure 5
#' Plot marginal coverage as function of the strength of label contamination,
#' increasing the number of labels
#' 

init_settings <- function(plot.optimistic = FALSE) {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+", "Label conditional+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)", "Adaptive+ (label-cond)")
    color.scale <<- cbPalette[c(1,3,4,8)]
    shape.scale <<- c(1,2,4,7)
    linetype.scale <<- c(1,1,1,1)
  }
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
}


load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num[1])
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
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.9,0.94), n_cal=1000, Method="Standard")
    
    {
      df3 = df2 = df[1:2,]
      df3$n_cal[1] = df2$n_cal[1] = min(df$n_cal)
      df3$n_cal[2] = df2$n_cal[2] = max(df$n_cal)
      df2$Mean[1] = 0.9
      df2$Mean[2] = 1.05
      df3$Mean[1] = 0.94
      df3$Mean[2] = 1.25
    }
    
    pp <- df %>%
      ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_point(data = df2, alpha = 0) +
      geom_point(data = df3, alpha = 0) +
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
      theme(text = element_text(size = 10),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            plot.margin = margin(5, -5, 1, -10)) 
    
    if(save_plots) {
      plot.file <- sprintf("figures/cifar10_%s_optimistic%s_%s_lc.pdf", plot.guarantee, plot.optimistic, plot.estimate)
      ggsave(file=plot.file, height=2.1, width=7, units="in")
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
        theme(text = element_text(size = 10),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/cifar10_%s_optimistic%s_%s_%d_lc.pdf",
                           plot.guarantee, plot.optimistic, plot.estimate, i)
      ggsave(file = plot.file, plot = pp, height = 3.5, width = 7, units = "in")
    }
  }
}


exp.num <- 101
plot.alpha <- 0.1
plot.K <- 10

# Figure 5
make_figure_101(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.estimate="rho-epsilon-point", plot.guarantee="marginal",
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)


### Experiment 201: BigEarthNet data ------------------------
#' Figure 6
#' Plot marginal coverage as function of the strength of label contamination, increasing the number of labels
#' 


init_settings <- function(plot.optimistic = FALSE) {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+", "Label conditional+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)", "Adaptive+ (label-cond)")
    color.scale <<- cbPalette[c(1,3,4,8)]
    shape.scale <<- c(1,2,4,7)
    linetype.scale <<- c(1,1,1,1)
  }
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
    group_by(data, K, n_cal, n_test, estimate,
             Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
  
  return(summary)
}


make_figure_201 <- function(exp.num, plot.alpha=0.1, plot.K, plot.estimate="rho-epsilon-point",
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
             Method %in% method.values, n_cal %in% c(500, 1500, 2500, 4500, 9500, 14500, 19500))%>%
      mutate(Method = factor(Method, method.values, method.labels),
             Mean = ifelse(Mean>1.8, NA, Mean))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.89,0.92), n_cal=1500, Method="Standard")
    
    pp <- df %>%
      ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.1) +
      facet_wrap(~Key, scales = "free_y", labeller = as_labeller(c("Coverage" = "Coverage", "Size" = "Size"))) +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      scale_x_continuous(trans='log10', limits=c(500,20000)) +
      scale_y_continuous(trans='log10') +
      xlab("Number of calibration samples") +
      ylab("") +
      theme_bw() +
      theme(text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12)) 
    
    if(save_plots) {
      plot.file <- sprintf("figures/bigearthnet_oracle_K%d_%s_optimistic%s_%s_lc.pdf",
                           plot.K, plot.guarantee, plot.optimistic, plot.estimate)
      ggsave(file=plot.file, height=3, width=7.5, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(Alpha==plot.alpha, K==plot.K, estimate==plot.estimate, Guarantee==plot.guarantee, Label=="marginal",
             Method %in% method.values, n_cal %in% c(500, 1500, 2500, 4500, 9500, 14500, 19000))%>%
      filter(n_cal >= 500)
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.89,0.92), n_cal=5000, Method="Standard")
    
    for (i in 1:length(method.values)) {
      current_methods <- method.values[1:i]
      current_labels <- method.labels[1:i]
      
      df_filtered <- df_filt %>%
        filter(Method %in% current_methods) %>%
        mutate(Method = factor(Method, levels = current_methods, labels = current_labels),
               Mean = ifelse(Mean>2.5, NA, Mean))
      
      {
        df3 = df2 = df[1:2,]
        df3$n_cal[1] = df2$n_cal[1] = min(df$n_cal)
        df3$n_cal[2] = df2$n_cal[2] = max(df$n_cal)
        df2$Mean[1] = 0.89
        df2$Mean[2] = 0.1
        df3$Mean[1] = 0.92
        df3$Mean[2] = 0.6
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
        scale_x_continuous(trans='log10', limits=c(500,20000)) +
        scale_y_continuous(trans='log10') +
        xlab("Number of calibration samples") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/bigearthnet_oracle_K%d_%s_optimistic%s_%s_lc_%d.pdf",
                           plot.K, plot.guarantee, plot.optimistic, plot.estimate, i)
      ggsave(file = plot.file, plot = pp, height = 3.5, width = 7, units = "in")
    }
  }
  
  
}


exp.num <- 201
plot.alpha <- 0.1
plot.K <- 6
plot.estimate <- "none"

# Figure 6
make_figure_201(exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
                plot.estimate="none", plot.guarantee="marginal",
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)


#' (Not shown in paper)
#' Plot label-conditional coverage as function of the strength of label contamination,
#' stratified by the number of labels
make_figure_202 <- function(exp.num, plot.alpha=0.1, plot.K, plot.estimate="none",
                            plot.guarantee="marginal",
                            plot.optimistic=FALSE, save_plots=FALSE, reload=TRUE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic=plot.optimistic)
  
  df <- summary %>%
    filter(Alpha==plot.alpha, K==plot.K, estimate==plot.estimate, Guarantee==plot.guarantee,
           Method %in% method.values, n_cal %in% c(500, 1500, 2500,4500, 9500), Label %in% label.values)

  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.98,0.92), n_cal=5000, Method="Standard")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Label = factor(Label, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    #        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
    facet_grid(Key~Label, scales="free") +
    #geom_point(data=df.range2, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/bigearthnet_oracle_K%d_%s_optimistic%s_%s_lcsize.pdf",
                         plot.K, plot.guarantee, plot.optimistic, plot.estimate)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 201
plot.alpha <- 0.1
plot.K <- 6
plot.estimate <- "none"

label.values <<- 0:5
label.labels <<- c("CWW", "Arable Land", "Agriculture", "Vegetation", "Urban Fabric", "Mixed")

## Figure A16
make_figure_202(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
                plot.estimate=plot.estimate,
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)




#' ---------------------------------------------------------------------------------------------------------------------
### Experiments 600: AP identification ------------------------
init_settings <- function() {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  
  # method.values <<- c("Clean sample", "benchmark", "EE", "IF", "LOF")
  # method.labels <<- c("Clean sample", "benchmark", "EE", "IF", "LOF")
  
  method.values <<- c("Clean sample", "SVC", "RFC", "EE", "IF", "optimal")
  method.labels <<- c("Clean sample", "SVC", "RFC", "EE", "IF", "Opt")
  
  color.scale <<- cbPalette[c(1,2,4,5,6,7)]
  shape.scale <<- c(1,0,3,4,5,6)
  linetype.scale <<- c(1,1,1,1,1,1)
}


#### Experiment 601: Impact of numerosity of the training set -----------------
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
    pivot_longer(c("size","accuracy","accuracy_tilde"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, pi_easy, contamination, flipy, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_601 <- function(exp.num, plot.data="syntheticAP", plot.K=4,
                            plot.pi_easy,
                            plot.contamination="uniform",
                            plot.flipy=0,
                            plot.epsilon=0.2,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==2, K==plot.K,
           pi_easy %in% plot.pi_easy,
           Method %in% method.values,
           contamination==plot.contamination,
           flipy==plot.flipy, epsilon==plot.epsilon)
  nominal_accuracy <- 1 - (plot.K-1)/plot.K*plot.flipy
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=nominal_accuracy)
  #df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Pi = sprintf("pi: %.2f", pi_easy)) %>%
    ggplot(aes(x=n_train1, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Pi, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    #geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the training set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_flipy%s_%s.png",
                         exp.num, plot.data, plot.K, plot.flipy, plot.contamination)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 601
plot.epsilon <- 0.1
plot.K <- 4
plot.contamination <- "uniform"
plot.pi_easy <- c(0.25, 0.5, 0.75, 1)

plot.flipy <- 0
plot.data <- "syntheticAP"
make_figure_601(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.pi_easy=plot.pi_easy,
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
                save_plots=TRUE, reload=TRUE)


#### Experiment 602: Impact of size of the set for T estimation -----------------
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
    pivot_longer(c("accuracy", "epsilon_res"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, pi_easy, contamination, flipy, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_602 <- function(exp.num, plot.data="syntheticAP", plot.K=4,
                            plot.pi_easy,
                            plot.contamination="uniform",
                            plot.flipy=0,
                            plot.epsilon=0.2,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==2, K==plot.K,
           pi_easy %in% plot.pi_easy,
           Method %in% method.values,
           contamination==plot.contamination,
           flipy==plot.flipy, epsilon==plot.epsilon)
  nominal_accuracy <- 1 - (plot.K-1)/plot.K*plot.flipy
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=nominal_accuracy)
  df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Pi = sprintf("pi: %.2f", pi_easy)) %>%
    ggplot(aes(x=n_train2, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Pi, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the anchor-selection set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_flipy%s_%s.png",
                         exp.num, plot.data, plot.K, plot.flipy, plot.contamination)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 602
plot.epsilon <- 0.1
plot.K <- 4
plot.contamination <- "uniform"
plot.pi_easy <- c(0.25, 0.5, 0.75, 1)

plot.flipy <- 0
plot.data <- "syntheticAP"
make_figure_602(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.pi_easy=plot.pi_easy,
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
                save_plots=TRUE, reload=TRUE)

#### Experiment 603: Impact of the contamination strength -----------------
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
    pivot_longer(c("size","accuracy"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, pi_easy, contamination, flipy, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_603 <- function(exp.num, plot.data="syntheticAP", plot.K=4,
                            plot.pi_easy=1,
                            plot.contamination="uniform",
                            plot.flipy=0,
                            plot.epsilon,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==2, K==plot.K,
           pi_easy==plot.pi_easy,
           Method %in% method.values,
           contamination==plot.contamination,
           flipy==plot.flipy, epsilon %in% plot.epsilon)
  nominal_accuracy <- 1 - (plot.K-1)/plot.K*plot.flipy
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=nominal_accuracy)
  #df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
    ggplot(aes(x=n_train1, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Epsilon, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    #geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the training set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_flipy%f_%s_pi%s.png",
                         exp.num, plot.data, plot.K, plot.flipy, plot.contamination, plot.pi_easy)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 603
plot.epsilon <- c(0, 0.05, 0.1, 0.2)
plot.K <- 4
plot.contamination <- "uniform"
plot.flipy <- 0
plot.data <- "syntheticAP"

plot.pi_easy <- 1
make_figure_603(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.pi_easy=plot.pi_easy,
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
                save_plots=TRUE, reload=TRUE)

plot.pi_easy <- 0.75
make_figure_603(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.pi_easy=plot.pi_easy,
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
                save_plots=TRUE, reload=TRUE)

#### Experiment 604: Impact of the contamination strength -----------------
# Performance indexes displayed as function of size of the set for T estimation
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
    pivot_longer(c("accuracy", "epsilon_res"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, pi_easy, contamination, flipy, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_604 <- function(exp.num, plot.data="syntheticAP", plot.K=4,
                            plot.pi_easy=1,
                            plot.contamination="uniform",
                            plot.flipy=0,
                            plot.epsilon,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==2, K==plot.K,
           pi_easy==plot.pi_easy,
           Method %in% method.values,
           contamination==plot.contamination,
           flipy==plot.flipy, epsilon %in% plot.epsilon)
  nominal_accuracy <- 1 - (plot.K-1)/plot.K*plot.flipy
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=nominal_accuracy)
  df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
    ggplot(aes(x=n_train2, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Epsilon, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the anchor-selection set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_flipy%s_%s_pi%s.png",
                         exp.num, plot.data, plot.K, plot.flipy, plot.contamination, plot.pi_easy)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 604
plot.epsilon <- c(0, 0.05, 0.1, 0.2)
plot.K <- 4
plot.contamination <- "uniform"
plot.flipy <- 0
plot.data <- "syntheticAP"

plot.pi_easy <- 1
make_figure_604(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.pi_easy=plot.pi_easy,
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
                save_plots=TRUE, reload=TRUE)

plot.pi_easy <- 0.75
make_figure_604(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.pi_easy=plot.pi_easy,
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
                save_plots=TRUE, reload=TRUE)

#### Experiment 605: Impact of the contamination strength -----------------
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
    pivot_longer(c("size","accuracy"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, pi_easy, center_scale, contamination, flipy, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_605 <- function(exp.num, plot.data="syntheticAP", plot.K=4,
                            plot.pi_easy=1,
                            plot.contamination="uniform",
                            plot.flipy=0,
                            plot.epsilon=0.1,
                            plot.center_scale,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==2, K==plot.K,
           pi_easy==plot.pi_easy,
           Method %in% method.values,
           contamination==plot.contamination,
           flipy==plot.flipy, epsilon==plot.epsilon,
           center_scale %in% plot.center_scale)
  nominal_accuracy <- 1 - (plot.K-1)/plot.K*plot.flipy
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=nominal_accuracy)
  #df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Center = sprintf("Sep: %.2f", center_scale)) %>%
    ggplot(aes(x=n_train1, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Center, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    #geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the training set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_flipy%f_%s_pi%s.png",
                         exp.num, plot.data, plot.K, plot.flipy, plot.contamination, plot.pi_easy)
    ggsave(file=plot.file, height=4.5, width=7.5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 605
plot.K <- 4
plot.contamination <- "uniform"
plot.flipy <- 0
plot.data <- "syntheticAP"
plot.pi_easy <- 1
plot.center_scale <- c(0.5, 0.75, 1)

plot.epsilon <- 0.1
make_figure_605(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.pi_easy=plot.pi_easy,
                plot.center_scale=plot.center_scale,
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
                save_plots=TRUE, reload=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
### Experiments 610: AP existence ------------------------
init_settings <- function() {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  
  # method.values <<- c("Split 5", "Split 10", "Split 20", "Boot 5", "Boot 10", "Boot 20")
  # method.labels <<- c("Split 5", "Split 10", "Split 20", "Boot 5", "Boot 10", "Boot 20")

  method.values <<- c("Split 5", "Boot 5")
  method.labels <<- c("Split", "Bootstrap")
  
  color.scale <<- cbPalette[c(1,2,4,5,6,7)]
  shape.scale <<- c(1,0,3,4,5,6)
  linetype.scale <<- c(1,1,1,1,1,1)
}


#### Experiment 611: Performances for different scenarios -----------------
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
    #pivot_longer(c("correct", "FP", "FN", "existence"), names_to="Key", values_to="Value") %>%
    #mutate(Key = recode(Key, correct="accuracy", FP="FPR", FN="FNR", existence="exist_rate")) %>%
    pivot_longer(c("existence"), names_to="Key", values_to="Value") %>%
    mutate(Key = recode(Key, existence="AP detection rate")) %>%
    group_by(data, scenario, contamination, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sqrt(mean(Value)*(1-mean(Value))/N), .groups = "drop")
    #summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N), .groups = "drop")
  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_611 <- function(exp.num, plot.data="syntheticAP",
                            plot.contamination="uniform",
                            plot.epsilon=0.1,
                            plot.n_train1=10000,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon,
           n_train1==plot.n_train1)
  
  df.exist_rate_line <- tibble(Key="AP detection rate",
                               Mean=c(1,0,0),
                               scenario=c("scenario1","scenario2","scenario3"))
  #df.exist_rate <- tibble(Key=c("AP detection rate","AP detection rate"), Mean=c(0,1), n_train2=1000, Method="Split 5")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Scenario = sprintf("%s", scenario)) %>%
    ggplot(aes(x=n_train2, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=pmax(0,Mean-SE), ymax=pmin(1,Mean+SE)), width = 0.1) +
    facet_grid(Key~Scenario, scales="free") +
    geom_hline(data=df.exist_rate_line, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the anchor-selection set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_eps%s_%s_nt1_%d.png",
                         exp.num, plot.data, plot.epsilon, plot.contamination, plot.n_train1)
    ggsave(file=plot.file, height=2.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 611
plot.contamination <- "uniform"
plot.data <- "syntheticAP"
plot.n_train1 <- 10000

plot.epsilon <- 0
make_figure_611(exp.num=exp.num, plot.data=plot.data,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon, plot.n_train1=plot.n_train1,
                save_plots=FALSE, reload=TRUE)

plot.epsilon <- 0.1
make_figure_611(exp.num=exp.num, plot.data=plot.data,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon, plot.n_train1=plot.n_train1,
                save_plots=TRUE, reload=TRUE)



### Experiments 620: T estimation with NN and EM ------------------------
init_settings <- function(sll_flag=FALSE) {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  
  
  # method.values <<- c("EM", "NN", "NN16", "NN SLL", "softmax")
  # method.labels <<- c("EM", "NN", "NN (slim)", "NN (sll)","softmax")
  if(sll_flag){
    # method.values <<- c("EM", "NN SLL", "NNw SLL", "NNuw SLL", "softmax")
    # method.labels <<- c("EM", "NNs", "NNs (weighted)", "NNs (upweighted)", "softmax")
    # color.scale <<- cbPalette[c(2,3,8,9,7)]
    # shape.scale <<- c(0,2,7,8,6)
    # linetype.scale <<- c(1,1,1,1,1)
    method.values <<- c("EM", "NN SLL", "NNw SLL", "softmax")
    method.labels <<- c("EM", "NNs", "NNs (weighted)", "softmax")
    color.scale <<- cbPalette[c(2,3,8,7)]
    shape.scale <<- c(0,2,7,6)
    linetype.scale <<- c(1,1,1,1)
  } else {
    # method.values <<- c("EM", "NN", "NNw", "NNuw", "softmax")
    # method.labels <<- c("EM", "NN", "NN (weighted)", "NN (upweighted)", "softmax")
    # color.scale <<- cbPalette[c(2,4,5,6,7)]
    # shape.scale <<- c(0,3,4,5,6)
    # linetype.scale <<- c(1,1,1,1,1)
    
    # method.values <<- c("EM", "NN", "NNw", "softmax")
    # method.labels <<- c("EM", "NN", "NN (weighted)", "softmax")
    # color.scale <<- cbPalette[c(2,4,5,7)]
    # shape.scale <<- c(0,3,4,6)
    # linetype.scale <<- c(1,1,1,1)
    
    method.values <<- c("EM", "NN", "softmax")
    method.labels <<- c("EM", "NN", "softmax")
    color.scale <<- cbPalette[c(2,4,7)]
    shape.scale <<- c(0,3,6)
    linetype.scale <<- c(1,1,1)
  }
  # method.values <<- c("EM", "NN", "NN SLL", "softmax")
  # method.labels <<- c("EM", "NN", "NN (sll)","softmax")
  
  # method.values <<- c("EM", "NN")
  # method.labels <<- c("EM", "NN")
  
  # color.scale <<- cbPalette[c(2,4,5,6,7)]
  # shape.scale <<- c(0,3,4,5,6)
  # linetype.scale <<- c(1,1,1,1,1)
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
    pivot_longer(c("epsilon_res", "frobenius_d", "accuracy"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, contamination, epsilon, n, n_clean, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

#### Experiment 621: Impact of size of clean data -----------------
#' Plot performance as function of the number of calibration samples,
#' increasing the number of clean data
#' The clean observations are "easy observations"
make_figure_621 <- function(exp.num, plot.data="synthetic6", plot.K=4,
                            plot.n_clean=100,
                            plot.contamination="uniform",
                            plot.epsilon=0.2,
                            plot.sll_flag=FALSE,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(sll_flag=plot.sll_flag)
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, K==plot.K,
           n_clean %in% plot.n_clean,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon)
  
  prova <- (df$Method=="softmax") & (df$Key %in% c("frobenius_d", "epsilon_res"))
  df$Mean[prova] <- NA
  
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=1)
  df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  df.nominal_res_dist <- tibble(Key="frobenius_d", Mean=0)
  df.range_accuracy <- tibble(Key=c("accuracy","accuracy"), Mean=c(0.5,1), n=1000, Method="EM")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(N_CLEAN = factor(sprintf("Size of clean data: %d", n_clean),
                            levels = sprintf("Size of clean data: %d", plot.n_clean),
                            labels = sprintf("Size of clean data: %d", plot.n_clean))) %>%
    # mutate(K_lab = factor(sprintf("%d classes", K), 
    #                       levels = sprintf("%d classes", c(10, 20, 50)), 
    #                       labels = c("10 classes", "20 classes", "50 classes"))) %>%
    ggplot(aes(x=n, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~N_CLEAN, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_res_dist, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range_accuracy, aes(x=n, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of training samples") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_%s.png",
                         exp.num, plot.data, plot.K, plot.contamination)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 621
plot.epsilon <- 0.1
plot.K <- 4
plot.contamination <- "uniform"
plot.n_clean <- c(100,500,1000,5000)
#plot.n_clean <- c(1000,5000, 10000, 20000)
plot.data <- "synthetic6"

make_figure_621(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.n_clean=plot.n_clean,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                plot.sll_flag=FALSE,
                save_plots=FALSE, reload=TRUE)

make_figure_621(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.n_clean=plot.n_clean,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                plot.sll_flag=TRUE,
                save_plots=FALSE, reload=TRUE)


#### Experiment 622: Impact of fraction of clean data -----------------
#' Plot performance indexes as function of the number of training samples,
#' increasing the fraction of clean data
#' The clean observations are "easy observations"
#' 

init_settings <- function() {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  
  method.values <<- c("EM", "NN", "NN SLL", "softmax")
  method.labels <<- c("EM", "NN", "NN (sll)","softmax")
  
  # method.values <<- c("EM", "NN")
  # method.labels <<- c("EM", "NN")
  
  color.scale <<- cbPalette[c(2,4,5,6,7)]
  shape.scale <<- c(0,3,4,5,6)
  linetype.scale <<- c(1,1,1,1,1)
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
    pivot_longer(c("epsilon_res", "frobenius_d", "accuracy"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, contamination, epsilon, n, pi_clean, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


make_figure_622 <- function(exp.num, plot.data="synthetic6", plot.K=4,
                            plot.pi_clean,
                            plot.contamination="uniform",
                            plot.epsilon=0.2,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, K==plot.K,
           pi_clean %in% plot.pi_clean,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon)
  
  prova <- (df$Method=="softmax") & (df$Key %in% c("frobenius_d", "epsilon_res"))
  df$Mean[prova] <- NA
  
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=1)
  df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  df.nominal_res_dist <- tibble(Key="frobenius_d", Mean=0)
  df.range_accuracy <- tibble(Key=c("accuracy","accuracy"), Mean=c(0.5,1), n=1000, Method="NN")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(PI_CLEAN = factor(sprintf("Frac. of clean data: %s", pi_clean),
                            levels = sprintf("Frac. of clean data: %s", plot.pi_clean),
                            labels = sprintf("Frac. of clean data: %s", plot.pi_clean))) %>%
    # mutate(K_lab = factor(sprintf("%d classes", K), 
    #                       levels = sprintf("%d classes", c(10, 20, 50)), 
    #                       labels = c("10 classes", "20 classes", "50 classes"))) %>%
    ggplot(aes(x=n, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~PI_CLEAN, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_res_dist, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range_accuracy, aes(x=n, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of training samples") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_%s.png",
                         exp.num, plot.data, plot.K, plot.contamination)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 622
plot.epsilon <- 0.1
plot.K <- 4
plot.contamination <- "uniform"
plot.pi_clean <- c(0.1,0.2,0.5,0.8)
plot.data <- "synthetic6"

make_figure_622(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.pi_clean=plot.pi_clean,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                save_plots=FALSE, reload=TRUE)


#### Experiment 623: Impact of contamination strength -----------------
#' Plot performance as function of the number of training samples,
#' increasing the contamination strength

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
    pivot_longer(c("epsilon_res", "frobenius_d", "accuracy"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, contamination, epsilon, n, n_clean, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


make_figure_623 <- function(exp.num, plot.data="synthetic6", plot.K=4,
                            plot.n_clean=100,
                            plot.contamination="uniform",
                            plot.epsilon=0.2,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, K==plot.K,
           epsilon %in% plot.epsilon,
           n_clean==plot.n_clean,
           Method %in% method.values,
           contamination==plot.contamination)
  
  prova <- (df$Method=="softmax") & (df$Key %in% c("frobenius_d", "epsilon_res"))
  df$Mean[prova] <- NA
  
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=1)
  df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  df.nominal_res_dist <- tibble(Key="frobenius_d", Mean=0)
  df.range_accuracy <- tibble(Key=c("accuracy","accuracy"), Mean=c(0.5,1), n=1000, Method="NN")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
    ggplot(aes(x=n, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Epsilon, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_res_dist, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range_accuracy, aes(x=n, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of training samples") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_%s.png",
                         exp.num, plot.data, plot.K, plot.contamination)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 623
plot.epsilon <- c(0, 0.05, 0.1, 0.2)
plot.K <- 4
plot.contamination <- "uniform"
plot.n_clean <- 100
plot.data <- "synthetic6"

make_figure_623(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.n_clean=plot.n_clean,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                save_plots=FALSE, reload=TRUE)


#### Experiment 624: Different data design -----------------
#' Plot performance as function of the number of training samples,
#' changing the data design
make_figure_624 <- function(exp.num, plot.data="synthetic6", plot.K=4,
                            plot.n_clean=100,
                            plot.contamination="uniform",
                            plot.epsilon=0.2,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data %in% plot.data, num_var==20, K==plot.K,
           epsilon==plot.epsilon,
           n_clean==plot.n_clean,
           Method %in% method.values,
           contamination==plot.contamination)
  
  prova <- (df$Method=="softmax") & (df$Key %in% c("frobenius_d", "epsilon_res"))
  df$Mean[prova] <- NA
  
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=1)
  df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  df.nominal_res_dist <- tibble(Key="frobenius_d", Mean=0)
  df.range_accuracy <- tibble(Key=c("accuracy","accuracy"), Mean=c(0.5,1), n=1000, Method="NN")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Data = sprintf("Data: %s", data)) %>%
    ggplot(aes(x=n, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Data, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_res_dist, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range_accuracy, aes(x=n, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of training samples") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_eps%s_ncl%d_K%d_%s.png",
                         exp.num, plot.epsilon, plot.n_clean, plot.K, plot.contamination)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 624
plot.epsilon <- 0.2
plot.K <- 4
plot.contamination <- "uniform"
plot.n_clean <- 100
plot.data <- c("synthetic1", "synthetic2", "synthetic3")

make_figure_624(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.n_clean=plot.n_clean,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                save_plots=TRUE, reload=TRUE)



### Experiments 700: Using the estimated T in the adaptive algorithm ------------------------
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
    group_by(data, num_var, K, signal, model_name, contamination, flipy, epsilon, estimate, n_train, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  
  # method.values <<- c("Standard", "Standard using AP", "Adaptive optimized+", "Adaptive optimized+ clean",
  #                     "Adaptive optimized+ AP D2L", "Adaptive optimized+ AP drop1", "Adaptive optimized+ AP drop05",
  #                     "Adaptive optimized+ AP param")
  # method.labels <<- c("Standard", "Standard (AP)", "Adaptive+", "Adaptive+ (clean)",
  #                     "Adaptive+ (AP D2L)", "Adaptive+ (AP drop 1%)", "Adaptive+ (AP drop 0.5%)",
  #                     "Adaptive+ (AP RRM)")
  # color.scale <<- cbPalette[c(1,2,3,4,5,6,7,8)]
  # shape.scale <<- c(1,0,2,3,4,5,6,7)
  # linetype.scale <<- c(1,1,1,1,1,1,1,1)
  
  method.values <<- c("Standard",
                      "Standard using AP",
                      "Adaptive optimized+",
                      "Adaptive optimized+ clean",
                      "Adaptive optimized+ AP")
                      #"Adaptive optimized+ AP param")
  method.labels <<- c("Standard",
                      "Standard (AP)",
                      "Adaptive+",
                      "Adaptive+ (clean)",
                      #"Adaptive+ (AP drop 1%)",
                      "Adaptive+ (AP)")
                      #"Adaptive+ (AP RRM)")
  color.scale <<- cbPalette[c(1,3,4,5,6)]
  shape.scale <<- c(1,2,3,4,5)
  linetype.scale <<- c(1,1,1,1,1)
}


#' ---------------------------------------------------------------------------------------------------------------------
#### Experiment 701: Impact of class separation ------------------------
#' Plot marginal coverage as function of the number of calibration samples, increasing the class separation
make_figure_701 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.flipy=0,
                          plot.epsilon,
                          plot.optimistic=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train==10000, K==plot.K, signal %in% plot.signal, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           flipy==plot.flipy, epsilon==plot.epsilon)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Signal = sprintf("Sep: %.2f", signal)) %>%
    #        mutate(Label = factor(Label, label.values, label.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Signal, scales="free") +
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
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          #legend.position = "bottom",
          #legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_flipy%s_%s_%s_optimistic%s.pdf",
                         exp.num, plot.data, 10000, plot.K, plot.flipy, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.2
plot.signal <- c(0.7,1.0,2.0)
plot.K <- 4
plot.contamination <- "uniform"

exp.num <- 701
plot.data <- "synthetic1_easy"
plot.flipy <- 0
make_figure_701(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
              save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

plot.flipy <- 0.01
make_figure_701(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
                save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

exp.num <- 702
plot.data <- "synthetic1"
plot.flipy <- 0
make_figure_701(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

plot.flipy <- 0.01
make_figure_701(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)



#### Experiment 703: Impact of size of sample used for AP identification ------------------------
#' Plot marginal coverage as function of the number of calibration samples, increasing the size of the sample used to identify
#' the anchor points and estimate T
#' 

make_figure_703 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                            plot.contamination="uniform",
                            plot.flipy=0,
                            plot.epsilon,
                            plot.optimistic=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train %in% plot.n_train, K==plot.K, signal==plot.signal, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           flipy==plot.flipy, epsilon==plot.epsilon)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(N_sample = factor(sprintf("N_T_estim: %d", n_train), 
                          levels = sprintf("N_T_estim: %d", plot.n_train), 
                          labels = c("N_T_estim: 5000", "N_T_estim: 10000",
                                     "N_T_estim: 50000", "N_T_estim: 100000"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~N_sample, scales="free") +
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
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          #legend.position = "bottom",
          #legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_flipy%s_%s_%s_optimistic%s.pdf",
                         exp.num, plot.data, 10000, plot.K, plot.flipy, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.2
plot.signal <- 1.0
plot.n_train <- c(5000, 10000, 50000, 100000)
plot.K <- 4
plot.contamination <- "uniform"

exp.num <- 703
plot.data <- "synthetic1"
# plot.flipy <- 0
# make_figure_703(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
#                 plot.contamination=plot.contamination,
#                 plot.flipy=plot.flipy, plot.epsilon=plot.epsilon, save_plots=FALSE, plot.optimistic=TRUE, reload=TRUE)

plot.flipy <- 0.01
make_figure_703(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
#### Experiment 704: Simulations with real Anchor Points ------------------------
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
    group_by(data, num_var, K, signal, contamination, flipy, epsilon, estimate, n_train1, n_train2, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  
  method.values <<- c("Standard",
                      "Standard using AP",
                      "Adaptive optimized+",
                      "Adaptive optimized+ clean",
                      "Adaptive optimized+ AP opt")
  #"Adaptive optimized+ AP param")
  method.labels <<- c("Standard",
                      "Standard (AP)",
                      "Adaptive+",
                      "Adaptive+ (clean)",
                      "Adaptive+ (AP)")
  #"Adaptive+ (AP RRM)")
  color.scale <<- cbPalette[c(1,3,4,5,6)]
  shape.scale <<- c(1,2,3,4,5)
  linetype.scale <<- c(1,1,1,1,1)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the class separation
make_figure_704a <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal",
                             plot.contamination="uniform",
                             plot.flipy=0,
                             plot.epsilon=0.1,
                             plot.n_train1,
                             plot.n_train2=10000,
                             plot.optimistic=TRUE,
                             save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data==plot.data, num_var==2, n_train2==plot.n_train2, K==plot.K, Guarantee==plot.guarantee,
           Label=="marginal", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           flipy==plot.flipy, epsilon==plot.epsilon, n_train1 %in% plot.n_train1)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Signal = factor(sprintf("N train: %d", n_train1), 
                           levels = sprintf("N train: %d", c(500, 1000, 5000, 10000)), 
                           labels = c("N train: 500", "N train: 1000", "N train: 5000", "N train: 10000"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Signal, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_nt2_%d_K%d_flipy%s_%s_%s_optimistic%s.pdf",
                         exp.num, plot.data, plot.n_train2, plot.K, plot.flipy, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 704
plot.data <- "syntheticAP"
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.K <- 4
plot.contamination <- "uniform"
plot.flipy <- 0


plot.n_train1 <- c(500,1000,5000,10000)
plot.n_train2 <- 500
make_figure_704a(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
                plot.contamination=plot.contamination,
                plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
                plot.n_train1=plot.n_train1, plot.n_train2=plot.n_train2,
                save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


make_figure_704b <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal",
                             plot.contamination="uniform",
                             plot.flipy=0,
                             plot.epsilon=0.1,
                             plot.n_train1=10000,
                             plot.n_train2,
                             plot.optimistic=TRUE,
                             save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data==plot.data, num_var==2, n_train1==plot.n_train1, K==plot.K, Guarantee==plot.guarantee,
           Label=="marginal", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           flipy==plot.flipy, epsilon==plot.epsilon, n_train2 %in% plot.n_train2)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Signal = factor(sprintf("N AP-sel. set: %d", n_train2), 
                          levels = sprintf("N AP-sel. set: %d", c(500, 1000, 5000, 10000)), 
                          labels = c("N AP-sel. set: 500", "N AP-sel. set: 1000", "N AP-sel. set: 5000", "N AP-sel. set: 10000"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Signal, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_nt1_%d_K%d_flipy%s_%s_%s_optimistic%s.pdf",
                         exp.num, plot.data, plot.n_train1, plot.K, plot.flipy, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.n_train1 <- 1000
plot.n_train2 <- c(500,1000,5000,10000)
make_figure_704b(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
                 plot.contamination=plot.contamination,
                 plot.flipy=plot.flipy, plot.epsilon=plot.epsilon,
                 plot.n_train1=plot.n_train1, plot.n_train2=plot.n_train2,
                 save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
### Experiments 800: AP identification in CIFAR-10 dataset ------------------------
init_settings <- function() {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  
  # method.values <<- c("Clean sample", "SVC", "ResNet18", "IF")
  # method.labels <<- c("Clean sample", "SVC", "ResNet18", "IF")
  
  # method.values <<- c("Clean sample", "SVC", "IF")
  # method.labels <<- c("Clean sample", "SVC", "IF")
  
  method.values <<- c("Clean sample", "SVC", "IF", "optimal")
  method.labels <<- c("Clean sample", "SVC", "IF", "Opt")
  
  color.scale <<- cbPalette[c(1,2,6,7)]
  shape.scale <<- c(1,0,5,6)
  linetype.scale <<- c(1,1,1,1)
}


#### Experiment 801: Impact of numerosity of the training set -----------------
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
    pivot_longer(c("size","accuracy","accuracy_tilde"), names_to = "Key", values_to = "Value") %>%
    group_by(data, K, contamination, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_801 <- function(exp.num,plot.data="cifar10",
                            plot.contamination="uniform",
                            plot.epsilon,
                            plot.n_train2=5000,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, Method %in% method.values,
           contamination==plot.contamination,
           epsilon%in%plot.epsilon,
           n_train2==plot.n_train2)

  df.nominal_accuracy <- tibble(Key="accuracy", Mean=1)
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
    ggplot(aes(x=n_train1, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Epsilon, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the training set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_%s_nt2_%d.png",
                         exp.num, plot.data, plot.contamination, plot.n_train2)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 801
plot.epsilon <- c(0, 0.05, 0.1, 0.15)
plot.contamination <- "uniform"
plot.n_train2 <- 1000

plot.data <- "cifar10"
make_figure_801(exp.num=exp.num, plot.data=plot.data,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                plot.n_train2=plot.n_train2,
                save_plots=TRUE, reload=TRUE)


#### Experiment 802: Impact of numerosity of the training set -----------------
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
    pivot_longer(c("accuracy", "epsilon_res"), names_to = "Key", values_to = "Value") %>%
    group_by(data, K, contamination, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_802 <- function(exp.num,plot.data="cifar10",
                            plot.contamination="uniform",
                            plot.epsilon,
                            plot.n_train1=5000,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, Method %in% method.values,
           contamination==plot.contamination,
           epsilon%in%plot.epsilon,
           n_train1==plot.n_train1)
  
  df.nominal_error <- tibble(Key="epsilon_res", Mean=0)
  #df.nominal_error2 <- tibble(Key="frobenius_d", Mean=0)
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
    ggplot(aes(x=n_train2, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Epsilon, scales="free") +
    geom_hline(data=df.nominal_error, aes(yintercept=Mean), linetype="dashed") +
    #geom_hline(data=df.nominal_error2, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the anchor-selection set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_%s_nt1_%d.png",
                         exp.num, plot.data, plot.contamination, plot.n_train1)
    ggsave(file=plot.file, height=3.5, width=7.5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 802
plot.epsilon <- c(0, 0.05, 0.1, 0.15)
plot.contamination <- "uniform"
plot.n_train1 <- 3000

plot.data <- "cifar10"
make_figure_802(exp.num=exp.num, plot.data=plot.data,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                plot.n_train1=plot.n_train1,
                save_plots=TRUE, reload=TRUE)


exp.num <- 1002
plot.epsilon <- c(0, 0.05, 0.1)
plot.contamination <- "uniform"
plot.n_train1 <- 4000

plot.data <- "bigearthnet"
make_figure_802(exp.num=exp.num, plot.data=plot.data,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                plot.n_train1=plot.n_train1,
                save_plots=TRUE, reload=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
### Experiments 810: AP existence in CIFAR dataset ------------------------
init_settings <- function() {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  
  # method.values <<- c("Split 5", "Split 10", "Split 20", "Boot 5", "Boot 10", "Boot 20")
  # method.labels <<- c("Split 5", "Split 10", "Split 20", "Boot 5", "Boot 10", "Boot 20")
  
  method.values <<- c("Split 5", "Boot 5")
  method.labels <<- c("Split", "Bootstrap")
  
  color.scale <<- cbPalette[c(1,2,4,5,6,7)]
  shape.scale <<- c(1,0,3,4,5,6)
  linetype.scale <<- c(1,1,1,1,1,1)
}


#### Experiment 812: Performances for different scenarios -----------------
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
    #pivot_longer(c("correct", "FP", "FN", "existence"), names_to="Key", values_to="Value") %>%
    #mutate(Key = recode(Key, correct="accuracy", FP="FPR", FN="FNR", existence="exist_rate")) %>%
    pivot_longer(c("existence"), names_to="Key", values_to="Value") %>%
    mutate(Key = recode(Key, existence="AP detection rate")) %>%
    group_by(data, contamination, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sqrt(mean(Value)*(1-mean(Value))/N), .groups = "drop")
    #summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N), .groups = "drop")
  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_812 <- function(exp.num, plot.data="syntheticAP",
                             plot.contamination="uniform",
                             plot.epsilon=0.1,
                             plot.n_train1=4000,
                             save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data,
           Method %in% method.values,
           contamination==plot.contamination,
           #epsilon==plot.epsilon,
           epsilon %in% plot.epsilon,
           n_train1==plot.n_train1)
  
  #df_exist_rate_line <- tibble(Key=c("exist_rate"), Mean=1, n_train2=1000, Method="Split 5")
  df.exist_rate <- tibble(Key=c("AP detection rate","AP detection rate"), Mean=c(0,1), n_train2=1000, Method="Split 5")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Epsilon = sprintf("%s", epsilon)) %>%
    ggplot(aes(x=n_train2, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=pmax(0,Mean-SE), ymax=pmin(1,Mean+SE)), width = 0.1) +
    facet_grid(Key~Epsilon, scales="free") +
    geom_point(data=df.exist_rate, aes(x=n_train2, y=Mean), alpha=0) +
    #geom_hline(data=df.exist_rate_line, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the anchor-selection set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_eps%s_%s_nt1_%d.png",
                         exp.num, plot.data, plot.epsilon, plot.contamination, plot.n_train1)
    ggsave(file=plot.file, height=2.5, width=5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 812
plot.epsilon <- c(0.15)
plot.contamination <- "uniform"
plot.data <- "cifar10"
plot.n_train1 <- 3000

make_figure_812(exp.num=exp.num, plot.data=plot.data,
                 plot.contamination=plot.contamination,
                 plot.epsilon=plot.epsilon, plot.n_train1=plot.n_train1,
                 save_plots=TRUE, reload=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
### Experiments 900: Noise-adaptive conformal in CIFAR-10 dataset ------------------------

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
    group_by(data, contamination, epsilon, n_train1, n_train2, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    filter(seed %in% (1:20)) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  
  method.values <<- c("Standard",
                      "Standard using AP",
                      "Adaptive optimized+",
                      "Adaptive optimized+ clean",
                      "Adaptive optimized+ AP SVC")
  #"Adaptive optimized+ AP param")
  method.labels <<- c("Standard",
                      "Standard (AP)",
                      "Adaptive+",
                      "Adaptive+ (clean)",
                      "Adaptive+ (AP)")
  #"Adaptive+ (AP RRM)")
  color.scale <<- cbPalette[c(1,3,4,5,6)]
  shape.scale <<- c(1,2,3,4,5)
  linetype.scale <<- c(1,1,1,1,1)
}


make_figure_901 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.guarantee="marginal",
                             plot.contamination="uniform",
                             plot.epsilon=0.1,
                             plot.nu=0,
                             plot.n_train1=1000,
                             plot.n_train2,
                             plot.optimistic=TRUE,
                             save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data==plot.data, n_train1==plot.n_train1, Guarantee==plot.guarantee,
           Label=="marginal", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, n_train2 %in% plot.n_train2)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.875,1), n_cal=1000, Method="Standard")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Signal = factor(sprintf("N AP-sel. set: %d", n_train2), 
                           levels = sprintf("N AP-sel. set: %d", c(1000, 2000, 3000)), 
                           labels = c("N AP-sel. set: 1000", "N AP-sel. set: 2000", "N AP-sel. set: 3000"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Signal, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of calibration samples") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_nt1_%d_eps%s_nu%s_%s_optimistic%s.pdf",
                         exp.num, plot.data, plot.n_train1, plot.epsilon, plot.nu, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 901
plot.data <- "cifar10"
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.contamination <- "uniform"

plot.n_train1 <- 1000
plot.n_train2 <- c(1000,2000,3000)
make_figure_901(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.guarantee="marginal",
                 plot.contamination=plot.contamination,
                 plot.epsilon=plot.epsilon,
                 plot.n_train1=plot.n_train1, plot.n_train2=plot.n_train2,
                 save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
### Experiments 1000: AP identification in BIGEARTHNET dataset ------------------------
init_settings <- function() {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  
  # method.values <<- c("Clean sample", "SVC", "ResNet18", "IF")
  # method.labels <<- c("Clean sample", "SVC", "ResNet18", "IF")
  
  # method.values <<- c("Clean sample", "SVC", "IF")
  # method.labels <<- c("Clean sample", "SVC", "IF")
  
  method.values <<- c("Clean sample", "SVC", "IF", "optimal")
  method.labels <<- c("Clean sample", "SVC", "IF", "Opt")
  
  color.scale <<- cbPalette[c(1,2,6,7)]
  shape.scale <<- c(1,0,5,6)
  linetype.scale <<- c(1,1,1,1)
}



#### Experiment 1002: Impact of numerosity of the anchor-selection set -----------------
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
    pivot_longer(c("accuracy", "epsilon_res"), names_to = "Key", values_to = "Value") %>%
    group_by(data, K, contamination, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_1002 <- function(exp.num,plot.data="cifar10",
                            plot.contamination="uniform",
                            plot.epsilon,
                            plot.n_train1=5000,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, Method %in% method.values,
           contamination==plot.contamination,
           epsilon%in%plot.epsilon,
           n_train1==plot.n_train1)
  
  df.nominal_error <- tibble(Key="epsilon_res", Mean=0)
  #df.nominal_error2 <- tibble(Key="frobenius_d", Mean=0)
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
    ggplot(aes(x=n_train2, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Epsilon, scales="free") +
    geom_hline(data=df.nominal_error, aes(yintercept=Mean), linetype="dashed") +
    #geom_hline(data=df.nominal_error2, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the anchor-selection set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_%s_nt1_%d.png",
                         exp.num, plot.data, plot.contamination, plot.n_train1)
    ggsave(file=plot.file, height=3.5, width=7.5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}


exp.num <- 1002
plot.epsilon <- c(0, 0.05, 0.1)
plot.contamination <- "uniform"
plot.n_train1 <- 4000

plot.data <- "bigearthnet"
make_figure_1002(exp.num=exp.num, plot.data=plot.data,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                plot.n_train1=plot.n_train1,
                save_plots=FALSE, reload=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
### Experiments 1010: AP existence in BEN dataset ------------------------
init_settings <- function() {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  
  # method.values <<- c("Split 5", "Split 10", "Split 20", "Boot 5", "Boot 10", "Boot 20")
  # method.labels <<- c("Split 5", "Split 10", "Split 20", "Boot 5", "Boot 10", "Boot 20")
  
  method.values <<- c("Split 5", "Boot 5")
  method.labels <<- c("Split", "Bootstrap")
  
  color.scale <<- cbPalette[c(1,2,4,5,6,7)]
  shape.scale <<- c(1,0,3,4,5,6)
  linetype.scale <<- c(1,1,1,1,1,1)
}


#### Experiment 1012: Performances for different scenarios -----------------
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
    #pivot_longer(c("correct", "FP", "FN", "existence"), names_to="Key", values_to="Value") %>%
    #mutate(Key = recode(Key, correct="accuracy", FP="FPR", FN="FNR", existence="exist_rate")) %>%
    pivot_longer(c("existence"), names_to="Key", values_to="Value") %>%
    mutate(Key = recode(Key, existence="AP detection rate")) %>%
    group_by(data, contamination, epsilon, n_train1, n_train2, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sqrt(mean(Value)*(1-mean(Value))/N), .groups = "drop")
    #summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N), .groups = "drop")
  
  return(summary)
}

#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_1012 <- function(exp.num, plot.data="syntheticAP",
                            plot.contamination="uniform",
                            plot.epsilon=0.1,
                            plot.n_train1=4000,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon %in% plot.epsilon,
           n_train1==plot.n_train1)
  
  #df_exist_rate_line <- tibble(Key=c("exist_rate"), Mean=1, n_train2=1000, Method="Split 5")
  df.exist_rate <- tibble(Key=c("AP detection rate","AP detection rate"), Mean=c(0,1), n_train2=1000, Method="Split 5")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Epsilon = sprintf("%s", epsilon)) %>%
    ggplot(aes(x=n_train2, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=pmax(0,Mean-SE), ymax=pmin(1,Mean+SE)), width = 0.1) +
    facet_grid(Key~Epsilon, scales="free") +
    geom_point(data=df.exist_rate, aes(x=n_train2, y=Mean), alpha=0) +
    #geom_hline(data=df.exist_rate_line, aes(yintercept=Mean), linetype="dashed") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of samples in the anchor-selection set") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_%s_nt1_%d.png",
                         exp.num, plot.data, plot.contamination, plot.n_train1)
    ggsave(file=plot.file, height=2.5, width=5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 1012
plot.epsilon <- c(0.15)
plot.contamination <- "uniform"
plot.data <- "bigearthnet"
plot.n_train1 <- 4000

make_figure_1012(exp.num=exp.num, plot.data=plot.data,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon, plot.n_train1=plot.n_train1,
                save_plots=TRUE, reload=TRUE)
