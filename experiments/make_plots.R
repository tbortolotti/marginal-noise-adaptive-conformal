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
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 11))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_nu%s_%s_%s_optimistic%s.pdf",
                           exp.num, plot.data, 10000, plot.K, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=4, width=6.5, units="in")
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
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=FALSE, slides=FALSE)

# Optimistic counterpart (not shown in paper)
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=FALSE, slides=FALSE)

# For slides
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=FALSE, slides=TRUE)
# Optimistic counterpart
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=FALSE, slides=TRUE)

plot.contamination <- "uniform"
## Figure A9
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=FALSE, slides=FALSE)


plot.contamination <- "block"
## Figure A10
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=FALSE, slides=FALSE)

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
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 11))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_K%d_eps%s_%s_%s_optimistic%s.pdf",
                           exp.num, 10000, plot.K, plot.epsilon, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=4, width=6.5, units="in")
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

exp.num <- 2

plot.alpha <- 0.1
plot.nu <- c(0, 0.25, 0.75, 1)
plot.epsilon <- 0.1
plot.K <- 4

## Figure 2
make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
# Optimistic counterpart (not shown in paper)
make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)


# For slides
make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=TRUE)
# Optimistic counterpart (not shown in paper)
make_figure_2(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=TRUE)


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
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 11))
    
    
    if(save_plots) {
      plot.file <- sprintf("figures/Aexp%d_synthetic1_ntrain%d_K%d_eps%s_%s_%s_optimistic%s.pdf",
                           exp.num, 10000, plot.K, plot.epsilon, plot.guarantee, plot.contamination, plot.optimistic)
      ggsave(file=plot.file, height=4, width=6.5, units="in")
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
make_figure_2A(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.ncal=plot.ncal, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)

#' ---------------------------------------------------------------------------------------------------------------------
### Experiment 3: Impact of number of classes, RRB -----------------------------
#' Figure 3 and Figure A14
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
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=6.5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}


exp.num <- 3
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.contamination <- "RRB"

## Figure 3
plot.nu <- 0.2
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# Optimistic counterpart (not shown in paper)
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

## Figure A14
plot.nu <- 0.8
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
## Figure 4
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
### Experiment 4: Impact of number of classes, Uniform and Block ------------------------
#' Figure A12 and Figure A13
#' Plot marginal coverage as function of the number of calibration samples, increasing the number of classes
#' 

exp.num <- 4
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.nu <- 0

## Figure A12: Uniform
plot.contamination <- "uniform"
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# Optimistic counterpart (not shown in paper)
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

## Figure A13: Block
plot.contamination <- "block"
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# Optimistic counterpart (not shown in paper)
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)


#' ---------------------------------------------------------------------------------------------------------------------
### Experiment 5: Impact of number of classes, RRB ------------------------
#' Figure A12 and Figure A13
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
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimistic%s.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=3.5, width=7, units="in")
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


## Figure 5
# Optimistic counterpart
make_figure_5(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

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
                          plot.optimistic=FALSE,
                          slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data=="synthetic4", num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, nu %in% plot.nu) %>%
    filter(n_cal >= 500)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  
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
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_K%d_eps%s_%s_%s_optimistic%s.pdf",
                         exp.num, 10000, plot.K, plot.epsilon, plot.guarantee, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=4, width=6.5, units="in")
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

## Figure 301
make_figure_301(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)


### Experiment 101: CIFAR-10H data ------------------------
#' Figure 6
#' Plot marginal coverage as function of the strength of label contamination,
#' increasing the number of labels
#' 

init_settings <- function(plot.optimistic = FALSE) {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)")
    color.scale <<- cbPalette[c(1,3,4)]
    shape.scale <<- c(1,2,4)
    linetype.scale <<- c(1,1,1)
  } else {
    method.values <<- c("Standard", "Adaptive simplified", "Adaptive optimized", "Asymptotic")
    method.labels <<- c("Standard", "Adaptive (simplified)", "Adaptive", "Adaptive (asymptotic)")
    color.scale <<- cbPalette[c(1,7,2,9)]
    shape.scale <<- c(1,3,0,5)
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
    
    {
      df3 = df2 = df[1:2,]
      df3$n_cal[1] = df2$n_cal[1] = min(df$n_cal)
      df3$n_cal[2] = df2$n_cal[2] = max(df$n_cal)
      df2$Mean[1] = 0.88
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
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 11)) 
    
    if(save_plots) {
      plot.file <- sprintf("figures/cifar10_%s_optimistic%s_%s.pdf", plot.guarantee, plot.optimistic, plot.estimate)
      ggsave(file=plot.file, height=3, width=6.5, units="in")
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

# Adaptive methods without the optimistic option (Not shown in paper)
make_figure_101(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.estimate="rho-epsilon-point", plot.guarantee="marginal",
                plot.optimistic=FALSE, save_plots=TRUE, reload=TRUE)
# Adaptive methods for marginal coverage with the optimistic option (Not shown in paper) 
make_figure_101(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.estimate="rho-epsilon-point", plot.guarantee="marginal",
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)

### Experiment 102: CIFAR-10H data ------------------------
#' Figure A15
#' Plot label-conditional coverage as function of the strength of label contamination,
#' stratified by the number of labels


init_settings <- function(plot.optimistic = FALSE) {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Label conditional+")
    method.labels <<- c("Standard", "Adaptive+ (label-cond)")
    color.scale <<- cbPalette[c(1,8)]
    shape.scale <<- c(1,7)
    linetype.scale <<- c(1,1)
  } else {
    method.values <<- c("Standard", "Label conditional")
    method.labels <<- c("Standard", "Adaptive (label-cond)")
    color.scale <<- cbPalette[c(1,10)]
    shape.scale <<- c(1,6)
    linetype.scale <<- c(1,1)
  }
  
}

make_figure_102 <- function(exp.num, plot.alpha=0.1, plot.K, plot.estimate="rho-epsilon-point",
                            plot.guarantee="marginal",
                            plot.optimistic=FALSE, save_plots=FALSE, reload=TRUE, fig.num=1) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic=plot.optimistic)
  
  df <- summary %>%
    filter(Alpha==plot.alpha, K==plot.K, estimate==plot.estimate, Guarantee==plot.guarantee,
           Method %in% method.values, n_cal %in% c(500,1500,4500,9500), Label %in% label.values) %>%
    filter(n_cal >= 500)
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
    plot.file <- sprintf("figures/cifar10_%s_optimistic%s_%s_lc_%d.pdf", plot.guarantee, plot.optimistic, plot.estimate,fig.num)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}


exp.num <- 102
plot.alpha <- 0.1
plot.K <- 10

label.values <<- 0:4
label.labels <<- paste("Label", 1:5, sep=" ")

## Figure A15 a)
make_figure_102(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
                plot.estimate="rho-epsilon-point",
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE, fig.num=1)

label.values <<- 5:9
label.labels <<- paste("Label", 6:10, sep=" ")

## Figure A15 b)
make_figure_102(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
                plot.estimate="rho-epsilon-point",
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE, fig.num=2)


# Experiment 101 and 102: CIFAR-10H ------------------------------------------------

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
  results1 <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  
  idir <- sprintf("results_hpc/exp%d", exp.num[2])
  ifile.list <- list.files(idir)
  results2 <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  idxs.out <- which(results2$Method=="Standard")
  
  results <- rbind(results1,results2[-idxs.out,])
  
  summary <- results %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, K, n_cal, n_test, epsilon_n_clean, epsilon_n_corr, estimate, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
  
  return(summary)
}

make_figure_103 <- function(exp.num, plot.alpha=0.1, plot.K, plot.estimate="rho-epsilon-point",
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
    
    {
      df3 = df2 = df[1:2,]
      df3$n_cal[1] = df2$n_cal[1] = min(df$n_cal)
      df3$n_cal[2] = df2$n_cal[2] = max(df$n_cal)
      df2$Mean[1] = 0.88
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
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 11)) 
    
    if(save_plots) {
      plot.file <- sprintf("figures/cifar10_%s_optimistic%s_%s_lcn.pdf", plot.guarantee, plot.optimistic, plot.estimate)
      ggsave(file=plot.file, height=3, width=6.5, units="in")
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
      
      plot.file <- sprintf("figures/slides/cifar10_%s_optimistic%s_%s_%d_lc.pdf",
                           plot.guarantee, plot.optimistic, plot.estimate, i)
      ggsave(file = plot.file, plot = pp, height = 3.5, width = 7, units = "in")
    }
  }
  
  
}


exp.num <- c(101,102)
plot.alpha <- 0.1
plot.K <- 10

# Figure 6
make_figure_103(exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.estimate="rho-epsilon-point", plot.guarantee="marginal",
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)


label.values <<- 0:4
label.labels <<- paste("Label", 1:5, sep=" ")

## Figure A15 a)
make_figure_102(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
                plot.estimate="rho-epsilon-point",
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE, fig.num=1)

label.values <<- 5:9
label.labels <<- paste("Label", 6:10, sep=" ")

## Figure A15 b)
make_figure_102(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
                plot.estimate="rho-epsilon-point",
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE, fig.num=2)




### Experiment 201: BigEarthNet data ------------------------
#' Figure 6
#' Plot marginal coverage as function of the strength of label contamination, increasing the number of labels
#' 

init_settings <- function(plot.optimistic = FALSE) {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  if(plot.optimistic) {
    #method.values <<- c("Standard", "Adaptive optimized+", "Adaptive simplified+", "Asymptotic+")
    #method.labels <<- c("Standard", "Adaptive-o+", "Adaptive-s+", "Asymptotic+")
    method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+")
    method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)")
    color.scale <<- cbPalette[c(1,3,4)]
    shape.scale <<- c(1,2,4)
    linetype.scale <<- c(1,1,1)
    # method.values <<- c("Standard", "Adaptive optimized+")
    # method.labels <<- c("Standard", "Adaptive+")
    # color.scale <<- cbPalette[c(1,3)]
    # shape.scale <<- c(1,2)
    # linetype.scale <<- c(1,1)
    
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
    group_by(data, K, n_cal, n_test, estimate,
             Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))
  
  return(summary)
}

make_figure_201 <- function(exp.num, plot.alpha=0.1, plot.K, plot.estimate="none",
                            plot.guarantee="marginal",
                            plot.optimistic=FALSE, save_plots=FALSE, reload=TRUE,
                            slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic=plot.optimistic)
  
  if(!slides){
    
    df <- summary %>%
      filter(Alpha==plot.alpha, K==plot.K, estimate==plot.estimate,
             Guarantee==plot.guarantee, Label=="marginal",
             Method %in% method.values, n_cal %in% c(500, 1500, 2500, 4500, 9500, 14500, 19500))  %>%
      mutate(Method = factor(Method, method.values, method.labels))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.89,0.92), n_cal=2500, Method="Standard")
    
    {
      df3 = df2 = df[1:2,]
      df3$n_cal[1] = df2$n_cal[1] = min(df$n_cal)
      df3$n_cal[2] = df2$n_cal[2] = max(df$n_cal)
      df2$Mean[1] = 0.89
      df2$Mean[2] = 1.3
      df3$Mean[1] = 0.92
      df3$Mean[2] = 1.4
    }
    
    
    pp <- df %>%
      ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
      geom_point() +
      geom_line() +
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.1) +
      geom_point(data = df2, alpha = 0) +
      geom_point(data = df3, alpha = 0) +
      facet_wrap(.~Key, scales="free") +
      geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
      geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
      scale_color_manual(values=color.scale) +
      scale_shape_manual(values=shape.scale) +
      scale_linetype_manual(values=linetype.scale) +
      scale_x_continuous(trans='log10', limits=c(500,20000)) +
      xlab("Number of calibration samples") +
      ylab("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 11)) 
    
    if(save_plots) {
      plot.file <- sprintf("figures/bigearthnet_oracle_K%d_%s_optimistic%s_%s.pdf",
                           plot.K, plot.guarantee, plot.optimistic, plot.estimate)
      #ggsave(file=plot.file, height=3.5, width=7, units="in")
      ggsave(file=plot.file, height=3, width=6.5, units="in")
      return(NULL)
    } else{
      return(pp)
    }
  } else {
    df_filt <- summary %>%
      filter(Alpha==plot.alpha, K==plot.K, estimate==plot.estimate,
             Guarantee==plot.guarantee, Label=="marginal",
             Method %in% method.values, n_cal %in% c(500, 1500, 2500,4500, 9500, 14500,19500))
    
    df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.89,0.93), n_cal=2500, Method="Standard")
    
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
        df2$Mean[2] = 1.3
        df3$Mean[1] = 0.94
        df3$Mean[2] = 1.4
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
        xlab("Number of calibration samples") +
        ylab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal")
      
      plot.file <- sprintf("figures/slides/bigearthnet_oracle_K%d_%s_optimistic%s_%s_%d.pdf",
                           plot.K, plot.guarantee, plot.optimistic, plot.estimate, i)
      ggsave(file = plot.file, plot = pp, height = 3.5, width = 7, units = "in")
    }
  }
}

exp.num <- 201
plot.alpha <- 0.1
plot.K <- 6
plot.estimate = "none"

# Adaptive methods without the optimistic option (not shown in paper)
make_figure_201(exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
                plot.estimate=plot.estimate, plot.guarantee="marginal",
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)

## Figure 7
make_figure_201(exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
                plot.estimate=plot.estimate, plot.guarantee="marginal",
                plot.optimistic=FALSE, save_plots=TRUE, reload=TRUE)

### Experiment 202: BigEarthNet data ------------------------
#' Figure A16
#' Plot label-conditional coverage as function of the strength of label contamination,
#' stratified by the number of labels


init_settings <- function(plot.optimistic = FALSE) {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2")
  if(plot.optimistic) {
    method.values <<- c("Standard", "Label conditional+")
    method.labels <<- c("Standard", "Adaptive+ (label-cond)")
    color.scale <<- cbPalette[c(1,8)]
    shape.scale <<- c(1,7)
    linetype.scale <<- c(1,1)
  } else {
    method.values <<- c("Standard", "Label conditional")
    method.labels <<- c("Standard","Adaptive (label-cond)")
    color.scale <<- cbPalette[c(1,10,1)]
    shape.scale <<- c(1,6)
    linetype.scale <<- c(1,1)
  }
  label.values <<- 0:5
  label.labels <<- c("CWW", "Arable Land", "Agriculture", "Vegetation", "Urban Fabric", "Mixed")
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
    summarise(Mean=mean(Value, na.rm=TRUE), N=sum(!is.na(Value)), SE=2*sd(Value, na.rm=TRUE)/sqrt(N))
  
  return(summary)
}

make_figure_202 <- function(exp.num, plot.alpha=0.1, plot.K, plot.estimate="none",
                            plot.guarantee="marginal",
                            plot.optimistic=FALSE, save_plots=FALSE, reload=TRUE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic=plot.optimistic)
  
  df <- summary %>%
    filter(Alpha==plot.alpha, K==plot.K, estimate==plot.estimate, Guarantee==plot.guarantee,
           Method %in% method.values, n_cal %in% c(2000,3000,5000,10000,15000), Label %in% label.values)
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  #df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=10000, Method="Standard")
  #df.range2 <- tibble(Key=c("Size","Size"), Mean=c(1,2), n_cal=1000, Method="")
  
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
    plot.file <- sprintf("figures/bigearthnet_oracle_K%d_%s_optimistic%s_%s_lc.pdf",
                         plot.K, plot.guarantee, plot.optimistic, plot.estimate)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 202
plot.alpha <- 0.1
plot.K <- 6
plot.estimate <- "none"

## Figure A16
make_figure_202(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K,
                plot.estimate=plot.estimate,
                plot.optimistic=TRUE, save_plots=TRUE, reload=TRUE)










