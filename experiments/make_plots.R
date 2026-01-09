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


### Experiments 600: Using the estimated T in the adaptive algorithm ------------------------
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
  label.values <<- c("10 classes", "20 classes", "50 classes")
  label.labels <<- c("10 classes", "20 classes", "50 classes")
  
  method.values <<- c("Standard", "Standard using AP", "Adaptive optimized+", "Adaptive optimized+ clean",
                      "Adaptive optimized+ AP D2L", "Adaptive optimized+ AP drop1", "Adaptive optimized+ AP drop05", "Adaptive optimized+ AP drop01",
                      "Adaptive optimized+ AP param")
  method.labels <<- c("Standard", "Standard (AP)", "Adaptive+", "Adaptive+ (clean)",
                      "Adaptive+ (AP D2L)", "Adaptive+ (AP drop 10%)", "Adaptive+ (AP drop 5%)", "Adaptive+ (AP drop 1%)",
                      "Adaptive+ (AP RRM)")
  color.scale <<- cbPalette[c(1,2,3,4,5,6,7,8,9)]
  shape.scale <<- c(1,0,2,3,4,5,6,7,8)
  linetype.scale <<- c(1,1,1,1,1,1,1,1,1)
}


#' ---------------------------------------------------------------------------------------------------------------------
#### Experiment 601: Impact of Label contamination strength ------------------------
#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_601 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
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
           nu==plot.nu, epsilon %in% plot.epsilon)
  
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
}

exp.num <- 601
plot.alpha <- 0.1
plot.nu <- 0
plot.epsilon <- c(0.05,0.1,0.2)
plot.K <- 4
plot.data <- "synthetic1"

plot.contamination <- "uniform"
make_figure_601(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=FALSE, plot.optimistic=TRUE, reload=TRUE)

plot.contamination <- "block"
make_figure_601(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


#' COMMENTO:
#' Cose da fare qui:
#' 1. Per large calibration sample sizes non dovrebbe esserci questa grande distanza tra la nominal
#'    e la empirical coverage. Usando il metodo che usa Matteo che tiene conto anche dell'incertezza legata
#'    alla stima di T si risolve questo problema? Cioè dovrei "tirare" su tutto.
#' 2. 

#' ---------------------------------------------------------------------------------------------------------------------
#### Experiment 602: Impact of Label contamination model ------------------------
#' Plot marginal coverage as function of the number of calibration samples, increasing the deviation of contamination
#' from a randomized response model
#' 

make_figure_602 <- function(exp.num, plot.alpha, plot.data, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon, plot.nu=0,
                          plot.optimistic=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  if(TRUE){
    df <- summary %>%
      filter(data==plot.data, num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
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
      plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_eps%s_%s_%s_optimistic%s.pdf",
                           exp.num, plot.data, 10000, plot.K, plot.epsilon, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=3.2, width=9, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  }
}

exp.num <- 602

plot.data <- "synthetic1"
plot.alpha <- 0.1
plot.nu <- c(0, 0.25, 0.75, 1)
plot.epsilon <- 0.1
plot.K <- 4

make_figure_602(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=FALSE, plot.optimistic=FALSE, reload=TRUE)
# Optimistic counterpart
make_figure_602(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=FALSE, plot.optimistic=TRUE, reload=TRUE)


exp.num <- 604
plot.data <- "synthetic2"
plot.alpha <- 0.1
plot.nu <- c(0, 0.25, 0.75, 1)
plot.epsilon <- 0.1
plot.K <- 4

make_figure_602(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
                plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=FALSE, plot.optimistic=FALSE, reload=TRUE)
# Optimistic counterpart
make_figure_602(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
                plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=FALSE, plot.optimistic=TRUE, reload=TRUE)


# #### Experiment 605: Does it work in problems with many classes? ------------------
# make_figure_605 <- function(exp.num=1, plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
#                           plot.contamination="uniform",
#                           plot.epsilon=0.1, plot.nu=0.2,
#                           plot.optimistic = FALSE) {
#   if(reload) {
#     summary <- load_data(exp.num)
#   }
#   
#   init_settings(plot.optimistic = plot.optimistic)
#   
#   df <- summary %>%
#     filter(data=="synthetic1", num_var==20, n_train==10000, Guarantee==plot.guarantee,
#            Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
#            Method %in% method.values,
#            contamination==plot.contamination,
#            epsilon==plot.epsilon, nu==plot.nu, n_cal>=500)
#   
#   df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
#   df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.85,1), n_cal=1000, Method="Standard")
#   
#   pp <- df %>%
#     # mutate(Method = factor(Method, method.values, method.labels),
#     #        Mean = ifelse(Key == "Size", log10(Mean), Mean),
#     #        SE = ifelse(Key == "Size", log10(Mean + SE) - log10(Mean), SE)) %>%
#     mutate(Method = factor(Method, method.values, method.labels)) %>%
#     # mutate(K_lab = factor(sprintf("%d classes", K), 
#     #                       levels = sprintf("%d classes", c(10, 20, 50)), 
#     #                       labels = c("10 classes", "20 classes", "50 classes"))) %>%
#     mutate(K_lab = factor(sprintf("%d classes", K),
#                           levels = sprintf("%d classes", c(10, 20)),
#                           labels = c("10 classes", "20 classes"))) %>%
#     ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
#     geom_point() +
#     geom_line() +
#     #        geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
#     #facet_grid(Key~K_lab, scales="free_y", labeller=custom_labeller)+
#     facet_wrap(Key~K_lab, scales="free_y")+
#     geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
#     geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
#     scale_color_manual(values=color.scale) +
#     scale_shape_manual(values=shape.scale) +
#     scale_linetype_manual(values=linetype.scale) +
#     #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
#     scale_x_continuous(trans='log10') +
#     scale_y_continuous(trans='log10')+
#     xlab("Number of calibration samples") +
#     ylab("") +
#     theme_bw() +
#     theme(text = element_text(size = 12),
#           axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#           legend.position = "bottom",
#           legend.direction = "horizontal",
#           legend.text = element_text(size = 12),
#           legend.title = element_text(size = 12),
#           plot.margin = margin(5, 5, 1, -10))
#   
#   if(save_plots) {
#     plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s.pdf",
#                          exp.num,
#                          10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
#     ggsave(file=plot.file, height=5, width=7.5, units="in")
#     #ggsave(file=plot.file, height=4, width=7.5, units="in")
#     return(NULL)
#   } else{
#     return(pp)
#   }
# }
# 
# 
# exp.num <- 605
# plot.alpha <- 0.1
# plot.epsilon <- 0.1
# plot.contamination <- "uniform"
# plot.nu <- 0
# 
# make_figure_605(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=FALSE, plot.optimistic=FALSE, reload=TRUE)
# make_figure_605(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#                 plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=FALSE, plot.optimistic=TRUE, reload=TRUE)
# 


### Experiments 700: Estimating the contamination process ------------------------
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
    pivot_longer(c("tv_d", "frobenius_d", "frob_inv_d"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, signal, model_name, contamination, epsilon, nu, gamma, n_train, n_cal, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

init_settings <- function() {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  method.values <<- c("Clean sample", "Clean sample (n_eq)", "Anchor points Patrini", "Anchor points empirical", "Anchor points empirical parametric")
  method.labels <<- c("Clean sample", "Clean sample (n_eq)", "AP (Patrini)", "AP (empirical)", "AP (empirical param)")
  # method.values <<- c("Clean sample")
  # method.labels <<- c("Clean sample")
  color.scale <<- cbPalette[c(1,2,3,4,5)]
  shape.scale <<- c(1,0,2,3,4)
  linetype.scale <<- c(1,1,1,1,1)
}


#' ---------------------------------------------------------------------------------------------------------------------
#### Experiment 701: Impact of Label contamination strength ------------------------
#' Plot marginal coverage as function of the number of calibration samples, increasing the contamination strength
make_figure_701 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4,
                            plot.contamination="uniform", plot.n_train=10000, plot.signal=1, plot.model_name="RFC",
                            plot.epsilon, plot.nu=0, plot.gamma=0.03,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train==plot.n_train, K==plot.K, signal==plot.signal,
           model_name==plot.model_name,
           Method %in% method.values,
           contamination==plot.contamination,
           nu==plot.nu, epsilon %in% plot.epsilon, gamma==plot.gamma)
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Epsilon = sprintf("Contam: %.2f", epsilon)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Epsilon, scales="free") +
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
    plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_nu%s_%s.pdf",
                         exp.num, plot.data, plot.n_train, plot.K, plot.nu, plot.contamination)
    ggsave(file=plot.file, height=7.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 701
plot.nu <- 0
plot.epsilon <- c(0.05,0.1,0.15,0.2)
plot.K <- 4
plot.data <- "synthetic1"
plot.contamination <- "uniform"
plot.n_train <- 10000
plot.signal <- 1
plot.model_name <- "RFC"
plot.gamma <- 0.03

make_figure_701(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
              plot.signal=plot.signal, plot.model_name=plot.model_name,
              plot.contamination=plot.contamination, plot.n_train=plot.n_train,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.gamma=plot.gamma,
              save_plots=TRUE, reload=TRUE)

# Qui dovrò pensare agli altri stimatori parametrici su altri esperimenti


#### Experiment 702: Impact of the threshold gamma  ------------------------
make_figure_702 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4,
                            plot.contamination="uniform", plot.n_train=10000, plot.signal=1, plot.model_name="RFC",
                            plot.epsilon=0.2, plot.nu=0, plot.gamma,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train==plot.n_train, K==plot.K, signal==plot.signal,
           model_name==plot.model_name,
           Method %in% method.values,
           contamination==plot.contamination,
           nu==plot.nu, epsilon==plot.epsilon, gamma %in% plot.gamma)
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Gamma = sprintf("Thresh: %.2f", gamma)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Gamma, scales="free") +
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
    plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_nu%s_%s.pdf",
                         exp.num, plot.data, plot.n_train, plot.K, plot.nu, plot.contamination)
    ggsave(file=plot.file, height=7.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 702
plot.nu <- 0
plot.epsilon <- 0.2
plot.K <- 4
plot.data <- "synthetic1"
plot.contamination <- "uniform"
plot.n_train <- 10000
plot.signal <- 1
plot.model_name <- "RFC"
plot.gamma <- c(0.01,0.05, 0.1, 0.2)

make_figure_702(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.signal=plot.signal, plot.model_name=plot.model_name,
                plot.contamination=plot.contamination, plot.n_train=plot.n_train,
                plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.gamma=plot.gamma,
                save_plots=TRUE, reload=TRUE)

#### Experiment 703: Impact of the class separation ------------------------
make_figure_703 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4,
                            plot.contamination="uniform", plot.n_train=10000, plot.signal, plot.model_name="RFC",
                            plot.epsilon=0.2, plot.nu=0, plot.gamma=0.03,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train==plot.n_train, K==plot.K,
           signal%in%plot.signal,
           model_name==plot.model_name,
           Method %in% method.values,
           contamination==plot.contamination,
           nu==plot.nu, epsilon==plot.epsilon, gamma==plot.gamma)
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(Signal = sprintf("Class sep: %.2f", signal)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~Signal, scales="free") +
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
    plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_nu%s_%s.pdf",
                         exp.num, plot.data, plot.n_train, plot.K, plot.nu, plot.contamination)
    ggsave(file=plot.file, height=7.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 703
plot.nu <- 0
plot.epsilon <- 0.2
plot.K <- 4
plot.data <- "synthetic1"
plot.contamination <- "uniform"
plot.n_train <- 10000
plot.signal <- c(0.1, 0.5, 1)
plot.model_name <- "RFC"
plot.gamma <- 0.03

make_figure_703(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.signal=plot.signal, plot.model_name=plot.model_name,
                plot.contamination=plot.contamination, plot.n_train=plot.n_train,
                plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.gamma=plot.gamma,
                save_plots=TRUE, reload=TRUE)


#### Experiment 704: Impact of the quality of the point predictor ------------------------
make_figure_704 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4,
                            plot.contamination="uniform", plot.n_train, plot.signal=1, plot.model_name="RFC",
                            plot.epsilon=0.2, plot.nu=0, plot.gamma=0.03,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train%in%plot.n_train, K==plot.K,
           signal==plot.signal,
           model_name==plot.model_name,
           Method %in% method.values,
           contamination==plot.contamination,
           nu==plot.nu, epsilon==plot.epsilon, gamma==plot.gamma)
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(N_train = factor(sprintf("N train: %d", n_train), 
                          levels = sprintf("N train: %d", c(500, 1000, 2000)), 
                          labels = c("N_train: 500", "N train: 1000", "N train: 2000"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~N_train, scales="free") +
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
    plot.file <- sprintf("figures/exp%d_%s_K%d_nu%s_%s.pdf",
                         exp.num, plot.data, plot.K, plot.nu, plot.contamination)
    ggsave(file=plot.file, height=7.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 704
plot.nu <- 0
plot.epsilon <- 0.2
plot.K <- 4
plot.data <- "synthetic1"
plot.contamination <- "uniform"
plot.n_train <- c(500, 1000, 2000)
plot.signal <- 1
plot.model_name <- "RFC"
plot.gamma <- 0.03

make_figure_704(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.signal=plot.signal, plot.model_name=plot.model_name,
                plot.contamination=plot.contamination, plot.n_train=plot.n_train,
                plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.gamma=plot.gamma,
                save_plots=TRUE, reload=TRUE)


#' Commento:
#' Probabilmente dovrei fare un unico plot stratificando per l'accuracy nell'individuare 
#' gli anchor points
#' Questi di fatto sono tutti plot che ci stanno dando la stessa informazione, ovvero che
#' quando perdo accuracy nel set di anchor points i metodi anchor points peggiorano.
#' L'altro commento che volevo fare è appunto che devo aggiungere altri stimatori, magari parametrici,
#' perché mi sembra strano il metodo Patrini et al. funzioni veramente così male.
#' Magari posso farmi venire in mente idee migliori, oltre al metodo ibrido già proposto?
#' 
#' Ci può stare che adesso io faccia qualche ragionamento sull'accuratezza.


#### Experiment 705: Analysis on the accuracy -----------------------------------------------------------------
#' Facciamo qualche prima visualizzazione...

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
    pivot_longer(c("accuracy", "accuracy_tilde", "frobenius_d", "tv_d"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, signal, model_name, contamination, epsilon, nu, gamma, n_train, n_cal, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}

init_settings <- function() {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  method.values <<- c("Anchor points Patrini", "Anchor points empirical", "Anchor points empirical parametric")
  method.labels <<- c("AP (Patrini)", "AP (empirical)", "AP (empirical param)")
  color.scale <<- cbPalette[c(2,3,4)]
  shape.scale <<- c(0,2,3)
  linetype.scale <<- c(1,1,1)
}

## Stratifying on n_cal
make_figure_705 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4,
                            plot.contamination="uniform", plot.n_train, plot.n_cal, plot.signal=1, plot.model_name="RFC",
                            plot.epsilon=0.2, plot.nu=0, plot.gamma,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train==plot.n_train, n_cal %in% plot.n_cal, K==plot.K,
           signal==plot.signal,
           model_name==plot.model_name,
           Method %in% method.values,
           contamination==plot.contamination,
           nu==plot.nu, epsilon==plot.epsilon, gamma %in% plot.gamma)
  
  df.range <- tibble(Key=c("accuracy","accuracy"), Mean=c(0.35,1), gamma=0.1, Method="AP (Patrini)")
  df.rangetilde <- tibble(Key=c("accuracy_tilde","accuracy_tilde"), Mean=c(0.35,1), gamma=0.1, Method="AP (Patrini)")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(N_cal = factor(sprintf("N cal: %d", n_cal), 
                          levels = sprintf("N cal: %d", c(500, 1000, 5000, 50000)), 
                          labels = c("N cal: 500", "N cal: 1000", "N cal: 5000", "N cal: 50000"))) %>%
    ggplot(aes(x=gamma, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_blank(data = df.range) +
    geom_blank(data = df.rangetilde) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~N_cal, scales="free") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10', limits=c(0.001,0.5)) +
    xlab("Threshold") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_nu%s_%s.pdf",
                         exp.num, plot.data, plot.K, plot.nu, plot.contamination)
    ggsave(file=plot.file, height=7.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 705
plot.nu <- 0
plot.epsilon <- 0.2
plot.K <- 4
plot.data <- "synthetic1"
plot.contamination <- "uniform"
plot.n_train <- 10000
plot.n_cal <- c(500, 1000, 5000, 50000)
plot.signal <- 1
plot.model_name <- "RFC"
plot.gamma <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)


make_figure_705(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.signal=plot.signal, plot.model_name=plot.model_name,
                plot.contamination=plot.contamination, plot.n_train=plot.n_train,
                plot.n_cal=plot.n_cal,
                plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.gamma=plot.gamma,
                save_plots=TRUE, reload=TRUE)

#### Experiment 706: Accuracy stratified on n_train -----------------------------------------------------------------

## Stratifying on n_train
make_figure_706 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4,
                            plot.contamination="uniform", plot.n_train, plot.n_cal, plot.signal=1, plot.model_name="RFC",
                            plot.epsilon=0.2, plot.nu=0, plot.gamma,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train %in% plot.n_train, n_cal== plot.n_cal, K==plot.K,
           signal==plot.signal,
           model_name==plot.model_name,
           Method %in% method.values,
           contamination==plot.contamination,
           nu==plot.nu, epsilon==plot.epsilon, gamma %in% plot.gamma)
  
  df.range <- tibble(Key=c("accuracy","accuracy"), Mean=c(1/plot.K,1), gamma=0.1, Method="AP (Patrini)")
  df.rangetilde <- tibble(Key=c("accuracy_tilde","accuracy_tilde"), Mean=c(1/plot.K,1), gamma=0.1, Method="AP (Patrini)")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(N_train = factor(sprintf("N train: %d", n_train), 
                          levels = sprintf("N train: %d", c(500, 1000, 10000)), 
                          labels = c("N train: 500", "N train: 1000", "N train: 10000"))) %>%
    ggplot(aes(x=gamma, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_blank(data = df.range) +
    geom_blank(data = df.rangetilde) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~N_train, scales="free") +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10', limits=c(0.001,0.5)) +
    xlab("Threshold") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_nu%s_%s.pdf",
                         exp.num, plot.data, plot.K, plot.nu, plot.contamination)
    ggsave(file=plot.file, height=7.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 706
plot.nu <- 0
plot.epsilon <- 0.2
plot.K <- 4
plot.data <- "synthetic1"
plot.contamination <- "uniform"
plot.n_train <- c(500, 1000,10000)
plot.n_cal <- 5000
plot.signal <- 1
plot.model_name <- "RFC"
plot.gamma <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)


make_figure_706(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.signal=plot.signal, plot.model_name=plot.model_name,
                plot.contamination=plot.contamination, plot.n_train=plot.n_train,
                plot.n_cal=plot.n_cal,
                plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.gamma=plot.gamma,
                save_plots=TRUE, reload=TRUE)



#### Experiment 707: Accuracy stratified on n_train -----------------------------------------------------------------
init_settings <- function() {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  method.values <<- c("Anchor points Patrini", "Anchor points empirical")
  method.labels <<- c("AP (Patrini)", "AP (empirical)")
  color.scale <<- cbPalette[c(2,3)]
  shape.scale <<- c(0,2)
  linetype.scale <<- c(1,1)
}

plot.contamination <- "block"
exp.num <- 707
make_figure_706(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.signal=plot.signal, plot.model_name=plot.model_name,
                plot.contamination=plot.contamination, plot.n_train=plot.n_train,
                plot.n_cal=plot.n_cal,
                plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.gamma=plot.gamma,
                save_plots=TRUE, reload=TRUE)

#' Questo esperimento ci dice che il minimo dell'errore commesso con il metodo AP empirical
#' corrisponde (giustamente) al punto in cui l'accuratezza del set anchor points comincia a crollare.
#' Questo ha senso perché per ottenere stime empiriche buone servono sia tanti punti
#' che questi punti siano degli anchor points accurati.
#' Come il numero di calibration samples aumenta, la necessità di arrivare fino al "limite" dell'accuratezza
#' si riduce perché a quel punto è possibile massimizzare l'accuratezza.
#' 
#' Qui secondo me stiamo toccando il discorso delle numerosità equivalenti.
#' Data una certa configurazione (data generation process, point classifier, altro?),
#' esiste una sorta di sample size equivalente a cui si raggiunge una sorta di plateu?
#' Come faccio ad accorgermene? Provo a guardare l'andamento dell'errore del solo metodo dei clean samples, per vedere se c'è un plateu.
#' No, questa storia del plateu è una cagata mi sa.
#' Allora qual è la questione qui? Ah forse posso provare calibration sampple sizes ancora maggiori e vedere che c'è una stabilizzazione su accuracy 1.


#' Ok, qui mi sto portando a casa due punti importanti:
#' 1. Patrini et al. non migliora con l'aumentare del calibration sample size,
#'    almeno con i point predictors che ho considerato io.
#'    Questo è sicuramente un problema che va sistemato.
#'    Il fatto che l'errore vada a 0 con il calibration sample size è una
#'    caratteristica imprescindibile.
#' 2. Patrini et al., quando ha un errore competitivo con il metodo empirico,
#'    è perché viene utilizzato su soglie molto piccole.
#' 3






