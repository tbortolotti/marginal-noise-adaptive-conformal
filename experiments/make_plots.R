options(width = 300)

library(tidyverse)
library(latex2exp)
library(RColorBrewer)


### Experiment 1: Impact of Label contamination strength -----------------------
#' Figure 1, Figure A9 and Figure A10
#' Plot marginal coverage as function of the number of calibration samples, increasing the strength
#' of the label contamination
#' 
load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
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
# make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
#              plot.contamination=plot.contamination,
#              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)

plot.contamination <- "uniform"
## Figure A16
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)


plot.contamination <- "block"
## Figure A17
make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
              plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)

# Optimistic counterpart (not shown in paper)
# make_figure_1(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K, plot.guarantee="marginal",
#              plot.contamination=plot.contamination,
#              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)

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

## Figure A18
make_figure_2A(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.ncal=plot.ncal, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE, slides=FALSE)
# Optimistic counterpart (not shown in paper)
# make_figure_2A(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
#               plot.ncal=plot.ncal, plot.epsilon=plot.epsilon, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE, slides=FALSE)


### Experiment 3: Impact of the number of classes ------------------------------------------------

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
    ggsave(file=plot.file, height=5, width=7.5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 3
plot.alpha <- 0.1
plot.nu <- 0.8

# Figure A19
plot.contamination <- "uniform"
plot.epsilon <- 0.2
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)

# Figure A20
plot.contamination <- "RRB"
plot.epsilon <- 0.1
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)

# Figure A21
plot.epsilon <- 0.2
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)

# Figure A22
plot.contamination <- "block"
plot.epsilon <- 0.1
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# # Optimistic counterpart (not shown in paper)
# make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

# Figure A23
plot.epsilon <- 0.2
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=FALSE, reload=TRUE)
# # Optimistic counterpart (not shown in paper)
# make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
#               plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

# The advantage of optimistic calibration
# Figure 3
plot.epsilon <- 0.2
plot.contamination <- "RRB"
plot.nu <- 0.8
make_figure_3(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)



### Experiment 4: The advantage of targeting marginal coverage 1 ------------------------------------------------
#' Contamination model: RRB with nu = 0.2 -------------------------------------
#' Plot marginal coverage as function of the number of calibration samples, increasing the number of classes
#' The separability of the classes increases as K increases, so that the average size of the
#' prediction sets remains stable across experiments with different K

load_data <- function(exp.num, plot.signal) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
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

make_figure_4 <- function(exp.num=1, plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
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


exp.num <- 4
plot.alpha <- 0.1
plot.epsilon <- 0.05
plot.contamination <- "RRB"
plot.nu <- 0.2
plot.signal <- 5

## Figure 4
make_figure_4(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, plot.signal=plot.signal, save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

### Experiment 201:The advantage of targeting marginal coverage 2 ------------------
#' RBB, Increase in the class-imbalance

load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
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


init_settings <- function() {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2", "#B22222")
  method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+", "Label conditional+")
  method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)", "Adaptive+ (label-cond)")
  color.scale <<- cbPalette[c(1,3,4,8)]
  shape.scale <<- c(1,2,4,7)
  linetype.scale <<- c(1,1,1,1)
}

make_figure_201 <- function(exp.num, plot.alpha, plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon, plot.nu=0,
                          imb.values,
                          plot.data="synthetic4") {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, n_train==10000, K==plot.K, signal==1, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, nu %in% plot.nu, imb %in% imb.values)
  
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
    plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_eps%s_nu%s_%s_%s_optimisticTRUE.pdf",
                         exp.num, plot.data, 10000, plot.K, plot.epsilon, plot.nu, plot.guarantee, plot.contamination)
    ggsave(file=plot.file, height=4, width=7.5, units="in")
    return(NULL)
  } else{
    return(pp)
  }
  
}


exp.num <- 201
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.contamination <- "RRB"
plot.nu <- 0.2
plot.K <- 4
plot.data <- "synthetic4"
imb.values <- c(0, 0.5, 1)

## Figure A24
make_figure_201(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, imb.values=imb.values, plot.data=plot.data, save_plots=TRUE, reload=TRUE)


### Experiment 202: Comparison with Clarkson (not shown in paper) -------------------
# We show what happens to the method of Clarkson et al. (2024) when the label frequencies
# are estimated from the calibration data, as the class imbalance increases

#' Class imbalance 
init_settings <- function() {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2", "#B22222")
  method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+", "Clarkson")
  method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)", "Clarkson et al.")
  color.scale <<- cbPalette[c(1,3,4,11)]
  shape.scale <<- c(1,2,4,8)
  linetype.scale <<- c(1,1,1,1)
}


load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
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

exp.num <- 202
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.contamination <- "RRB"
plot.nu <- 0.2
plot.K <- 10
plot.data <- "synthetic4"
imb.values <- c(0, 0.5, 1)

## Not shown in paper
#make_figure_201(exp.num=exp.num, plot.alpha=plot.alpha, plot.K=plot.K, plot.guarantee="marginal", plot.contamination="RRB",
#                plot.epsilon=plot.epsilon, plot.nu=plot.nu, imb.values=imb.values, plot.data=plot.data, save_plots=TRUE, reload=TRUE)



### Experiment 203: Comparison with Clarkson ------------------------------------
#' Increasing number of classes 

init_settings <- function() {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2", "#B22222")
  method.values <<- c("Standard", "Adaptive optimized+", "Asymptotic+", "Clarkson")
  method.labels <<- c("Standard", "Adaptive+", "Adaptive+ (asymptotic)", "Clarkson et al.")
  color.scale <<- cbPalette[c(1,3,4,11)]
  shape.scale <<- c(1,2,4,8)
  linetype.scale <<- c(1,1,1,1)
}

load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
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

make_figure_203 <- function(exp.num=1, plot.alpha=0.1, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon=0.1, plot.nu=0.2) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data=="synthetic1", num_var==20, n_train==10000, signal==1, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, nu==plot.nu)
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
          #legend.position = "bottom",
          #legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 5, 1, -10))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_synthetic1_ntrain%d_eps%f_nu%s_%s_%s_optimisticTRUE.pdf",
                         exp.num,
                         10000, plot.epsilon, plot.nu, plot.guarantee, plot.contamination)
    ggsave(file=plot.file, height=3.2, width=8, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values <<- c("4 classes", "8 classes", "16 classes")
label.labels <<- c("4 classes", "8 classes", "16 classes")
exp.num <- 203
plot.alpha <- 0.1
plot.epsilon <- 0.1
plot.nu <- 0.2

# Figure A7
plot.contamination <- "uniform"
make_figure_203(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
              plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, reload=TRUE)

# Figure A9
plot.contamination <- "RRB"
make_figure_203(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, reload=TRUE)

# Figure A11
plot.contamination <- "block"
make_figure_203(exp.num=exp.num, plot.alpha=plot.alpha, plot.guarantee="marginal", plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon, plot.nu=plot.nu, save_plots=TRUE, reload=TRUE)



### Experiment 204: Comparison with Clarkson 3 ----------------------------------
#' Increasing epsilon
make_figure_204 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                          plot.contamination="uniform",
                          plot.epsilon, plot.nu=0,
                          slides=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
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
    plot.file <- sprintf("figures/exp%d_%s_ntrain%d_K%d_nu%s_%s_%s_optimisticTRUE.pdf",
                         exp.num, plot.data, 10000, plot.K, plot.nu, plot.guarantee, plot.contamination)
    ggsave(file=plot.file, height=3.2, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
  
}

exp.num <- 204
plot.alpha <- 0.1
plot.nu <- 0.2
plot.epsilon <- c(0.05,0.1,0.2,0.4)
plot.K <- 4
plot.data <- "synthetic1"

# Figure A8
plot.contamination <- "uniform"
make_figure_204(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
              plot.guarantee="marginal",
              plot.contamination=plot.contamination, plot.epsilon=plot.epsilon, plot.nu=plot.nu,
              save_plots=TRUE, reload=TRUE, slides=FALSE)

# Figure A10
plot.contamination <- "RRB"
make_figure_204(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.guarantee="marginal",
                plot.contamination=plot.contamination, plot.epsilon=plot.epsilon, plot.nu=plot.nu,
                save_plots=TRUE, reload=TRUE, slides=FALSE)

# Figure A12
plot.contamination <- "block"
make_figure_204(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.guarantee="marginal",
                plot.contamination=plot.contamination, plot.epsilon=plot.epsilon, plot.nu=plot.nu,
                save_plots=TRUE, reload=TRUE, slides=FALSE)


### Experiments 300: T estimation with NN ------------------------
init_settings <- function(sll_flag=FALSE) {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  
  method.values <<- c("NN SLL alt",
                      "NN alt",
                      "softmax")
  method.labels <<- c("NNs",
                      "NN",
                      "softmax")
  color.scale <<- cbPalette[c(3,5,8,9,6)]
  shape.scale <<- c(3,7,4,6,5)
  linetype.scale <<- c(1,1,1,1,1)
}

#### Experiment 301: Impact of size of clean data -----------------
#' Plot performance as function of the number of calibration samples,
#' increasing the number of clean data
#' The clean observations are "easy observations"
#' 

load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE) 
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))    
  summary <- results %>%
    pivot_longer(c("epsilon_res", "accuracy"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, contamination, epsilon, n, n_clean, n_noisy, randflag, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


make_figure_301 <- function(exp.num, plot.data="synthetic6", plot.K=4,
                            plot.n_clean=100,
                            plot.randflag,
                            plot.contamination="uniform",
                            plot.epsilon=0.2,
                            save_plots=FALSE, reload=FALSE,sll_flag=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(sll_flag=sll_flag)
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, K==plot.K,
           n_clean %in% plot.n_clean,
           randflag==plot.randflag,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon)
  
  prova <- (df$Method=="softmax") & (df$Key %in% c("frobenius_d", "epsilon_res"))
  df$Mean[prova] <- NA
  
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=1)
  df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  #df.nominal_res_dist <- tibble(Key="frobenius_d", Mean=0)
  df.range_accuracy <- tibble(Key=c("accuracy","accuracy"), Mean=c(0.5,1), n_noisy=1000, Method="NN")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(N_CLEAN = factor(sprintf("Size of clean data: %d", n_clean),
                            levels = sprintf("Size of clean data: %d", plot.n_clean),
                            labels = sprintf("Size of clean data: %d", plot.n_clean))) %>%
    # mutate(K_lab = factor(sprintf("%d classes", K), 
    #                       levels = sprintf("%d classes", c(10, 20, 50)), 
    #                       labels = c("10 classes", "20 classes", "50 classes"))) %>%
    ggplot(aes(x=n_noisy, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~N_CLEAN, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    #geom_hline(data=df.nominal_res_dist, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range_accuracy, aes(x=n_noisy, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of noisy training samples") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_%s_%s.png",
                         exp.num, plot.data, plot.K, plot.contamination, plot.randflag)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 301
plot.epsilon <- 0.2
plot.K <- 4
plot.contamination <- "uniform"
plot.n_clean <- c(100,500,1000)
plot.data <- "synthetic6"

# Figure A28
make_figure_301(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.n_clean=plot.n_clean,
                plot.randflag=FALSE,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                save_plots=TRUE, reload=TRUE)


#### Experiment 302: Impact of fraction of clean data -----------------
#' Plot performance indexes as function of the number of training samples,
#' increasing the fraction of clean data
#' The clean observations are "easy observations"
#' 

load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE) 
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))    
  summary <- results %>%
    pivot_longer(c("epsilon_res", "accuracy"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, contamination, epsilon, n, pi_clean, randflag, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


make_figure_302 <- function(exp.num, plot.data="synthetic6", plot.K=4,
                            plot.pi_clean,
                            plot.rand_flag,
                            plot.contamination="uniform",
                            plot.epsilon=0.2,
                            save_plots=FALSE, reload=FALSE, sll_flag=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.sll_flag)
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, K==plot.K,
           pi_clean %in% plot.pi_clean,
           randflag==plot.rand_flag,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon)
  
  prova <- (df$Method=="softmax") & (df$Key %in% c("frobenius_d", "epsilon_res"))
  df$Mean[prova] <- NA
  
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=1)
  df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  #df.nominal_res_dist <- tibble(Key="frobenius_d", Mean=0)
  df.range_accuracy <- tibble(Key=c("accuracy","accuracy"), Mean=c(0.5,1), n=1000, Method="EM")
  
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
    #geom_hline(data=df.nominal_res_dist, aes(yintercept=Mean), linetype="dashed") +
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

exp.num <- 302
plot.epsilon <- 0.2
plot.K <- 4
plot.contamination <- "uniform"
plot.pi_clean <- c(0.1,0.3,0.5)
plot.data <- "synthetic6"

## Not shown in paper
#make_figure_302(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
#                plot.pi_clean=plot.pi_clean,
#                plot.rand_flag=FALSE,
#                plot.contamination=plot.contamination,
#                plot.epsilon=plot.epsilon,
#                save_plots=TRUE, reload=TRUE)

#### Experiment 303: Impact of contamination strength -----------------
#' Plot performance as function of the number of training samples,
#' increasing the contamination strength

load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE) 
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))    
  summary <- results %>%
    pivot_longer(c("epsilon_res"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, contamination, epsilon, n, n_clean, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


make_figure_303 <- function(exp.num, plot.data="synthetic6", plot.K=4,
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
    ggsave(file=plot.file, height=2.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 303
plot.epsilon <- c(0, 0.05, 0.1, 0.2)
plot.K <- 4
plot.contamination <- "uniform"
plot.n_clean <- 100
plot.data <- "synthetic6"

## Not shown in paper
#make_figure_303(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
#                plot.n_clean=plot.n_clean,
#                plot.contamination=plot.contamination,
#                plot.epsilon=plot.epsilon,
#                save_plots=TRUE, reload=TRUE)


#### Experiment 304: Different data design -----------------
#' Plot performance as function of the number of training samples,
#' changing the data design
#' 

load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE) 
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))    
  summary <- results %>%
    pivot_longer(c("epsilon_res"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, contamination, epsilon, n, n_clean, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}



make_figure_304 <- function(exp.num, plot.data="synthetic6", plot.K=4,
                            plot.n_clean=100,
                            plot.contamination="uniform",
                            plot.epsilon=0.2,
                            save_plots=FALSE, reload=FALSE, sll_flag=TRUE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(sll_flag)
  
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
    #geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
    #geom_hline(data=df.nominal_res_dist, aes(yintercept=Mean), linetype="dashed") +
    #geom_point(data=df.range_accuracy, aes(x=n, y=Mean), alpha=0) +
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
    ggsave(file=plot.file, height=2.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 304
plot.epsilon <- 0.2
plot.K <- 4
plot.contamination <- "uniform"
plot.n_clean <- 500
plot.data <- c("synthetic1", "synthetic2", "synthetic3")

# Figure A29
make_figure_304(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                plot.n_clean=plot.n_clean,
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                save_plots=TRUE, reload=TRUE)



#### Experiment 305: Different contamination model -----------------
#' Plot performance as function of the number of training samples,
#' changing the contamination model

init_settings <- function(sll_flag=FALSE) {
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#8A2BE2", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#F0E442")
  method.values <<- c("NN alt gen", "NN SLL alt gen", "softmax")
  method.labels <<- c("NN", "NNs", "softmax")
  color.scale <<- cbPalette[c(9,6,8)]
  shape.scale <<- c(9,3,7)
  linetype.scale <<- c(1,1,1)
  
}

load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE) 
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))    
  summary <- results %>%
    pivot_longer(c("frobenius_d", "accuracy"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, contamination, epsilon, n, n_clean, randflag, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


make_figure_305 <- function(exp.num, plot.data="synthetic6", plot.K=4,
                             plot.n_clean=500,
                             plot.contamination,
                             plot.epsilon=0.2,
                             save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20, K==plot.K,
           epsilon==plot.epsilon,
           n_clean==plot.n_clean,
           Method %in% method.values,
           contamination %in% plot.contamination)
  
  prova <- (df$Method=="softmax") & (df$Key %in% c("frobenius_d", "epsilon_res"))
  df$Mean[prova] <- NA
  
  df.nominal_accuracy <- tibble(Key="accuracy", Mean=1)
  #df.nominal_residual <- tibble(Key="epsilon_res", Mean=0)
  df.nominal_res_dist <- tibble(Key="frobenius_d", Mean=0)
  df.range_accuracy <- tibble(Key=c("accuracy","accuracy"), Mean=c(0.5,1), n=1000, Method="NN")
  
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    #mutate(CONT = sprintf("Cont.: %s", contamination)) %>%
    mutate(CONT = factor(sprintf("Cont: %s", contamination),
                         levels = sprintf("Cont: %s", plot.contamination),
                         labels = c("Cont: block", "Cont: two-level", "Cont: near-diag"))) %>%
    ggplot(aes(x=n, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~CONT, scales="free") +
    geom_hline(data=df.nominal_accuracy, aes(yintercept=Mean), linetype="dashed") +
    #geom_hline(data=df.nominal_residual, aes(yintercept=Mean), linetype="dashed") +
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
    plot.file <- sprintf("figures/exp%d_%s_eps%s_ncl%s_K%d.png",
                         exp.num, plot.data, plot.epsilon, plot.n_clean, plot.K)
    ggsave(file=plot.file, height=4.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 305
plot.epsilon <- 0.2
plot.K <- 4
plot.contamination <- c( "block", "RRB", "mild")
plot.n_clean <- 500
plot.data <- "synthetic6"

# Figure A30
make_figure_305(exp.num=exp.num, plot.data=plot.data, plot.K=plot.K,
                 plot.n_clean=plot.n_clean,
                 plot.contamination=plot.contamination,
                 plot.epsilon=plot.epsilon,
                 save_plots=TRUE, reload=TRUE)

### Experiments 400: Using the estimated T in the adaptive algorithm ------------------------
load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE) 
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))    
  summary <- results %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, num_var, K, model_name, contamination, epsilon, estimate, n_train, n_clean, pi_clean, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2","#648767")
  
  method.values <<- c("Standard",
                      "Standard using clean",
                      "Adaptive optimized+",
                      "Adaptive optimized+ clean",
                      "Adaptive optimized+ NN",
                      "__spacer__",
                      "Standard (clean) line")
  #"Adaptive optimized+ AP param")
  method.labels <<- c("Standard",
                      "Standard (clean)",
                      "Adaptive+",
                      "Adaptive+ (c/n)",
                      "Adaptive+ (NN)",
                      "",
                      "Standard (clean, simple)")
  #"Adaptive+ (AP RRM)")
  color.scale <<- cbPalette[c(1,2,3,5,NA,10)]
  shape.scale <<- c(1,5,2,4,NA,NA)
  linetype.scale <<- c(1,1,1,1,0,4)
  xlab_ <<- "Number of noisy calibration samples"
}


#' ---------------------------------------------------------------------------------------------------------------------
#### Experiment 401: Impact of the number of clean samples ------------------------
#' Plot marginal coverage as function of the number of calibration samples, increasing the size of the
#' clean dataset

make_figure_401 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                            plot.contamination="uniform",
                            plot.n_train, plot.n_clean, plot.pi_clean,
                            plot.epsilon,
                            plot.exp_easy=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.exp_easy)
  
  df <- summary %>%
    filter(data==plot.data, num_var==20,
           n_train==plot.n_train, n_clean %in% plot.n_clean,
           #pi_clean==plot.pi_clean,
           K==plot.K, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon)
  
  df.clean.values <- df %>%
    filter(Method == "Standard using clean") %>%
    group_by(Key, n_clean) %>%                        # <-- stratify by n_clean
    summarise(Mean = mean(Mean), .groups = "drop")
  
  df.clean.legend <- df %>%
    group_by(Key, n_clean) %>%                        # <-- group by both
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal), .groups = "drop") %>%
    left_join(df.clean.values, by = c("Key", "n_clean")) %>%   # <-- join on both
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "Standard (clean) line", SE = 0)
  
  df.spacer.legend <- df %>%
    group_by(Key, n_clean) %>%                        # <-- group by both
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal), .groups = "drop") %>%
    left_join(df.clean.values, by = c("Key", "n_clean")) %>%   # <-- join on both
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "__spacer__", SE = 0)

  
  # Aggiunge la riga fittizia al df principale (senza Standard using clean)
  df.plot <- df %>%
    filter(Method != "Standard using clean") %>%
    bind_rows(df.clean.legend) %>%
    bind_rows(df.spacer.legend)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.875,1), n_cal=1000, Method="Adaptive+")
  pp <- df.plot %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(N_CLEAN = factor(sprintf("N clean: %d", n_clean), 
                             levels = sprintf("N clean: %d", plot.n_clean), 
                             labels = c("N clean: 100", "N clean: 500",
                                        "N clean: 1000"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~N_CLEAN, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    xlab(xlab_) +
    ylab("") +
    guides(                                                                # <--
      color    = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1))),
      shape    = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1))),
      linetype = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1)))
    ) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          #legend.position = "bottom",
          #legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_nt%d_%s_%s_expeasy%s.pdf",
                         exp.num, plot.data, plot.K, plot.n_train, plot.guarantee, plot.contamination, plot.exp_easy)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.2
plot.n_train <- 5000
plot.n_clean <- c(100, 500, 1000)
plot.pi_clean <- 0
plot.K <- 4
plot.contamination <- "uniform"
exp.num <- 401
plot.data <- "synthetic6"

# Figure 5
make_figure_401(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.guarantee="marginal",
                plot.contamination=plot.contamination,
                plot.n_train=plot.n_train, plot.n_clean=plot.n_clean,
                plot.pi_clean=plot.pi_clean,
                plot.epsilon=plot.epsilon, save_plots=TRUE, reload=TRUE)

#### Experiment 402: Impact of the fraction of clean samples ------------------------
#' Plot marginal coverage as function of the number of calibration samples, increasing the fraction of the
#' clean dataset

make_figure_402 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                            plot.contamination="uniform",
                            plot.n_train, plot.pi_clean,
                            plot.epsilon) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20,
           n_train==plot.n_train,
           pi_clean %in% plot.pi_clean,
           K==plot.K, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Adaptive+")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(PI_CLEAN = factor(sprintf("Frac. of clean data: %s", pi_clean), 
                             levels = sprintf("Frac. of clean data: %s", plot.pi_clean), 
                             labels = c("Frac. of clean data: 0.05", "Frac. of clean data: 0.1",
                                        "Frac. of clean data: 0.2"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~PI_CLEAN, scales="free") +
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
    plot.file <- sprintf("figures/exp%d_%s_K%d_nt%d_%s_%s.pdf",
                         exp.num, plot.data, plot.K, plot.n_train, plot.guarantee, plot.contamination)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.2
plot.n_train <- 10000
plot.pi_clean <- c(0.05, 0.1, 0.2)
plot.K <- 4
plot.contamination <- "uniform"
exp.num <- 402
plot.data <- "synthetic6"

## Not shown in paper
#make_figure_402(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
#                plot.guarantee="marginal",
#                plot.contamination=plot.contamination,
#                plot.n_train=plot.n_train, plot.pi_clean=plot.pi_clean,
#                plot.epsilon=plot.epsilon, save_plots=TRUE, reload=TRUE)

#### Experiment 403: Impact of the number of training samples ------------------------
#' Plot marginal coverage as function of the number of calibration samples,
#' increasing the size of the training set

make_figure_403 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                            plot.contamination="uniform",
                            plot.n_train, plot.n_clean, plot.pi_clean,
                            plot.epsilon, plot.exp_easy=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.exp_easy)
  
  if(plot.exp_easy){
    df <- summary %>%
      filter(data==plot.data, num_var==20,
             n_train==5000, n_clean==plot.n_clean,
             #pi_clean==plot.pi_clean,
             K==plot.K, Guarantee==plot.guarantee,
             Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
             Method %in% method.values,
             contamination==plot.contamination,
             epsilon==plot.epsilon)
    df.range <- tibble(Key=c("Coverage","Coverage", "Size", "Size"), Mean=c(0.8,0.925,1,2), n_cal=1000, Method="Standard")
  } else {
    df <- summary %>%
      filter(data==plot.data, num_var==20,
             n_train %in% plot.n_train, n_clean==plot.n_clean,
             #pi_clean==plot.pi_clean,
             K==plot.K, Guarantee==plot.guarantee,
             Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
             Method %in% method.values,
             contamination==plot.contamination,
             epsilon==plot.epsilon)
    df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.875,1), n_cal=1000, Method="Adaptive+")
  }
  
  df.clean.values <- df %>%
    filter(Method == "Standard using clean") %>%
    group_by(Key, n_train) %>%
    summarise(Mean = mean(Mean), .groups = "drop")
  
  df.clean.legend <- df %>%
    group_by(Key, n_train) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal), .groups = "drop") %>%
    left_join(df.clean.values, by = c("Key", "n_train")) %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "Standard (clean) line", SE = 0)
  
  df.spacer.legend <- df %>%
    group_by(Key, n_train) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal), .groups = "drop") %>%
    left_join(df.clean.values, by = c("Key", "n_train")) %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "__spacer__", SE = 0)
  
  # Aggiunge la riga fittizia al df principale (senza Standard using clean)
  df.plot <- df %>%
    filter(Method != "Standard using clean") %>%
    bind_rows(df.clean.legend) %>%
    bind_rows(df.spacer.legend)

  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  pp <- df.plot %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(N_TRAIN = factor(sprintf("N train: %d", n_train), 
                            levels = sprintf("N train: %d", plot.n_train), 
                            labels = c("N train: 1000", "N train: 5000",
                                       "N train: 10000"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1)
  
  if(plot.exp_easy){
    pp <- pp +
      facet_wrap(.~Key, scales="free")
  }else{
    pp <- pp +
      facet_grid(Key~N_TRAIN, scales="free")
  }
  pp <- pp +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab(xlab_) +
    ylab("") +
    guides(                                                                # <--
      color    = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1))),
      shape    = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1))),
      linetype = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1)))
    ) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          #legend.position = "bottom",
          #legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_ncl%d_%s_%s_expeasy%s.pdf",
                         exp.num, plot.data, plot.K, plot.n_clean, plot.guarantee, plot.contamination, plot.exp_easy)
    
    if(plot.exp_easy){
      ggsave(file=plot.file, height=2.5, width=7, units="in")
    }else{
      ggsave(file=plot.file, height=4, width=9, units="in")
    }
    
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.2
plot.n_train <- c(1000, 5000, 10000)
plot.pi_clean <- 0
plot.K <- 4
plot.contamination <- "uniform"
exp.num <- 403
plot.data <- "synthetic6"

# Figure A25
plot.n_clean <- 500
make_figure_403(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.guarantee="marginal",
                plot.contamination=plot.contamination,
                plot.n_train=plot.n_train, plot.n_clean=plot.n_clean,
                plot.pi_clean=plot.pi_clean,
                plot.epsilon=plot.epsilon, save_plots=TRUE, reload=TRUE)

# Figure A26
plot.n_clean <- 100
make_figure_403(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.guarantee="marginal",
                plot.contamination=plot.contamination,
                plot.n_train=plot.n_train, plot.n_clean=plot.n_clean,
                plot.pi_clean=plot.pi_clean,
                plot.epsilon=plot.epsilon, save_plots=TRUE, reload=TRUE)

#### Experiment 404: Impact of the number of training samples ------------------------
#' Plot marginal coverage as function of the number of calibration samples,
#' increasing the size of the training sample

make_figure_404 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                            plot.contamination="uniform",
                            plot.n_train, plot.pi_clean,
                            plot.epsilon) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings()
  
  df <- summary %>%
    filter(data==plot.data, num_var==20,
           n_train %in% plot.n_train,
           pi_clean==plot.pi_clean,
           K==plot.K, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.8,1), n_cal=1000, Method="Standard")
  pp <- df %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(N_TRAIN = factor(sprintf("N train: %d", n_train), 
                            levels = sprintf("N train: %d", plot.n_train), 
                            labels = c("N train: 1000", "N train: 5000",
                                       "N train: 10000"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~N_TRAIN, scales="free") +
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
    plot.file <- sprintf("figures/exp%d_%s_K%d_picl_%s_%s_%s.pdf",
                         exp.num, plot.data, plot.K, plot.pi_clean, plot.guarantee, plot.contamination)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.2
plot.n_train <- c(1000, 5000, 10000)
plot.pi_clean <- 0.1
plot.K <- 4
plot.contamination <- "uniform"
exp.num <- 404
plot.data <- "synthetic6"

## Not shown in paper
#make_figure_404(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
#                plot.guarantee="marginal",
#                plot.contamination=plot.contamination,
#                plot.n_train=plot.n_train, plot.pi_clean=plot.pi_clean,
#                plot.epsilon=plot.epsilon, save_plots=FALSE, reload=TRUE)
#

#' ---------------------------------------------------------------------------------------------------------------------
#### Experiment 405: Impact of different data distribution ------------------------
#' Plot marginal coverage as function of the number of calibration samples,
#' changing the data distribution


make_figure_405 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                            plot.contamination="uniform",
                            plot.n_train, plot.n_clean, plot.pi_clean,
                            plot.epsilon,
                            plot.exp_easy=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.exp_easy)
  
  df <- summary %>%
    filter(data%in%plot.data,
           n_train==plot.n_train, n_clean==plot.n_clean,
           K==plot.K, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon)
  
  df.clean.values <- df %>%
    filter(Method == "Standard using clean") %>%
    group_by(Key, data) %>%
    summarise(Mean = mean(Mean), .groups = "drop")
  
  df.clean.legend <- df %>%
    group_by(Key, data) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal), .groups = "drop") %>%
    left_join(df.clean.values, by = c("Key", "data")) %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "Standard (clean) line", SE = 0)
  
  df.spacer.legend <- df %>%
    group_by(Key, data) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal), .groups = "drop") %>%
    left_join(df.clean.values, by = c("Key", "data")) %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "__spacer__", SE = 0)
  
  # Aggiunge la riga fittizia al df principale (senza Standard using clean)
  df.plot <- df %>%
    filter(Method != "Standard using clean") %>%
    bind_rows(df.clean.legend) %>%
    bind_rows(df.spacer.legend)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.875,1), n_cal=1000, Method="Adaptive+")
  pp <- df.plot %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(DATA = sprintf("Data: %s", data)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~DATA, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    xlab(xlab_) +
    ylab("") +
    guides(                                                                # <--
      color    = guide_legend(override.aes = list(alpha = c(1,1,1,1,1,0,1))),
      shape    = guide_legend(override.aes = list(alpha = c(1,1,1,1,1,0,1))),
      linetype = guide_legend(override.aes = list(alpha = c(1,1,1,1,1,0,1)))
    ) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          #legend.position = "bottom",
          #legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_K%d_nt%d_ncl%d_%s_%s_expeasy%s.pdf",
                         exp.num, plot.K, plot.n_train, plot.n_clean, plot.guarantee, plot.contamination, plot.exp_easy)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.2
plot.n_train <- 10000
plot.n_clean <- 500
plot.pi_clean <- 0
plot.K <- 4
plot.contamination <- "uniform"
exp.num <- 405
plot.data <- c("synthetic1","synthetic2","synthetic3")

## Not shown in paper
#make_figure_405(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
#                plot.guarantee="marginal",
#                plot.contamination=plot.contamination,
#                plot.n_train=plot.n_train, plot.n_clean=plot.n_clean,
#                plot.pi_clean=plot.pi_clean,
#                plot.epsilon=plot.epsilon, save_plots=TRUE, reload=TRUE)

#' ---------------------------------------------------------------------------------------------------------------------
#### Experiment 406: Impact of the contamination process ------------------------
#' Plot marginal coverage as function of the number of calibration samples,
#' changing the contamination process
#' 
#' 


make_figure_406 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.K=4, plot.guarantee="marginal", save_plots=FALSE, reload=FALSE,
                            plot.contamination,
                            plot.n_train, plot.n_clean, plot.pi_clean,
                            plot.epsilon,
                            plot.exp_easy=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.exp_easy)
  
  if(plot.data=="synthetic1_easy"){
    num_var = 2
  } else {
    num_var = 2
  }
  
  df <- summary %>%
    filter(data==plot.data, num_var==num_var,
           n_train==plot.n_train, n_clean==plot.n_clean,
           K==plot.K, Guarantee==plot.guarantee,
           Label=="marginal", model_name=="RFC", Alpha==plot.alpha,
           Method %in% method.values,
           contamination %in% plot.contamination,
           epsilon==plot.epsilon)
  
  df.clean.values <- df %>%
    filter(Method == "Standard using clean") %>%
    group_by(Key, contamination) %>%
    summarise(Mean = mean(Mean), .groups = "drop")
  
  df.clean.legend <- df %>%
    group_by(Key, contamination) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal), .groups = "drop") %>%
    left_join(df.clean.values, by = c("Key", "contamination")) %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "Standard (clean) line", SE = 0)
  
  df.spacer.legend <- df %>%
    group_by(Key, contamination) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal), .groups = "drop") %>%
    left_join(df.clean.values, by = c("Key", "contamination")) %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "__spacer__", SE = 0)
  
  # Aggiunge la riga fittizia al df principale (senza Standard using clean)
  df.plot <- df %>%
    filter(Method != "Standard using clean") %>%
    bind_rows(df.clean.legend) %>%
    bind_rows(df.spacer.legend)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.875,1), n_cal=1000, Method="Adaptive+")
  pp <- df.plot %>%
    mutate(Method = factor(Method, method.values, method.labels),
           Mean = ifelse((Key == "Coverage") & (Mean<0.8), NA, Mean)) %>%
    mutate(CONT = factor(sprintf("Cont: %s", contamination),
                         levels = sprintf("Cont: %s", plot.contamination),
                         labels = c("Cont: block", "Cont: two-level", "Cont: near-diag"))) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_grid(Key~CONT, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    #        scale_x_continuous(trans='log10', breaks=c(1000,2000,5000,10000,20000)) +
    scale_x_continuous(trans='log10') +
    xlab(xlab_) +
    ylab("") +
    guides(                                                                # <--
      color    = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1))),
      shape    = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1))),
      linetype = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1)))
    ) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          #legend.position = "bottom",
          #legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_K%d_nt%d_ncl%d_%s_expeasy%s.pdf",
                         exp.num, plot.data, plot.K, plot.n_train, plot.n_clean, plot.guarantee, plot.exp_easy)
    ggsave(file=plot.file, height=4, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

plot.alpha <- 0.1
plot.epsilon <- 0.2
plot.n_train <- 10000
plot.n_clean <- 500
plot.pi_clean <- 0
plot.K <- 4
plot.contamination <- c("block", "RRB", "mild")
exp.num <- 406
plot.data <- "synthetic6"

# Figure A27
make_figure_406(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.K=plot.K,
                plot.guarantee="marginal",
                plot.contamination=plot.contamination,
                plot.n_train=plot.n_train, plot.n_clean=plot.n_clean,
                plot.pi_clean=plot.pi_clean,
                plot.epsilon=plot.epsilon, save_plots=TRUE, reload=TRUE)

#' ---------------------------------------------------------------------------------------------------------------------
### Experiment 503: Noise-adaptive conformal in CIFAR-10 dataset ------------------------
load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE) 
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))    
  summary <- results %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, contamination, epsilon, n_train, n_clean, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2","#648767")
  
  method.values <<- c("Standard",
                      "Standard using clean",
                      "Adaptive optimized+ NN",
                      "Asymptotic+ NN",
                      "Label conditional+",
                      "__spacer__",
                      "Standard (clean) line")
  
  method.labels <<- c("Standard",
                      "Standard (clean)",
                      "Adaptive+",
                      "Adaptive+ (asymptotic)",
                      "Adaptive+ (label-cond)",
                      "",
                      "Standard (clean, simple)")
  color.scale <<- cbPalette[c(1,3,4,7,NA,10)]
  shape.scale <<- c(1,2,3,5,NA,NA)
  linetype.scale <<- c(1,1,1,1,0,4)
}

make_figure_503 <- function(exp.num, plot.alpha, plot.data="synthetic1", plot.guarantee="marginal",
                            plot.contamination="uniform",
                            plot.epsilon=0.1,
                            plot.nu=0,
                            plot.n_train=1000,
                            plot.n_clean,
                            plot.optimistic=TRUE,
                            save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data==plot.data, n_train==plot.n_train, Guarantee==plot.guarantee,
           Label=="marginal", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, n_clean %in% plot.n_clean)
  
  df.clean.values <- df %>%
    filter(Method=="Standard using clean") %>%
    group_by(Key) %>%
    summarise(mean_values=mean(Mean))
  df.clean.coverage <- as.numeric(df.clean.values[1,2])
  df.clean.size <- as.numeric(df.clean.values[2,2])
  df.clean <- tibble(Key=c("Coverage","Size"), Mean=c(df.clean.coverage,df.clean.size))
  
  df.clean.legend <- df %>%
    group_by(Key) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal)) %>%
    left_join(df.clean, by = "Key") %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "Standard (clean) line", SE = 0)
  
  df.spacer.legend <- df %>%
    group_by(Key) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal)) %>%
    left_join(df.clean, by = "Key") %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "__spacer__", SE = 0)
  
  # Aggiunge la riga fittizia al df principale (senza Standard using clean)
  df.plot <- df %>%
    filter(Method != "Standard using clean") %>%
    bind_rows(df.clean.legend) %>%
    bind_rows(df.spacer.legend)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.89,0.92), n_cal=1000, Method="Standard")
  
  pp <- df.plot %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_wrap(~Key, scales = "free_y", labeller = as_labeller(c("Coverage" = "Coverage", "Size" = "Size"))) +
    #facet_grid(Key~N_CLEAN, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of noisy calibration samples") +
    ylab("") +
    guides(
      color    = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1))),
      shape    = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1))),
      linetype = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1)))
    ) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_nt%d_ncl%d_eps%s_%s_optimistic%s.pdf",
                         exp.num, plot.data, plot.n_train, plot.n_clean, plot.epsilon, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=2.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 503
plot.data <- "cifar10"
plot.alpha <- 0.1
plot.epsilon <- 0.051
plot.contamination <- "real"
plot.n_train <- 2000
plot.n_clean <- c(500)

# Figure A31
make_figure_503(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.guarantee="marginal",
                plot.contamination=plot.contamination,
                plot.epsilon=plot.epsilon,
                plot.n_train=plot.n_train, plot.n_clean=plot.n_clean,
                save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

#' ---------------------------------------------------------------------------------------------------------------------
### Experiment 601: Noise-adaptive conformal in BigEarthNet dataset ------------------------
load_data <- function(exp.num) {
  idir <- sprintf("results_hpc/exp%d", exp.num)
  ifile.list <- list.files(idir, recursive = FALSE) 
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))    
  summary <- results %>%
    pivot_longer(c("Coverage", "Size"), names_to = "Key", values_to = "Value") %>%
    group_by(data, contamination, epsilon, n_train, n_clean, n_cal, Guarantee, Alpha, Label, Method, Key) %>%
    summarise(Mean=mean(Value), N=n(), SE=2*sd(Value)/sqrt(N))  
  return(summary)
}


init_settings <- function(plot.optimistic = FALSE) {
  df.dummy <<- tibble(key="Coverage", value=0.95)
  df.dummy2 <<- tibble(key="Coverage", value=0.5)
  cbPalette <<- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#20B2AA", "#8A2BE2","#648767")
  
  method.values <<- c("Standard",
                      "Standard using clean",
                      "Adaptive+ NN",
                      "Asymptotic+ NN",
                      "Label conditional+",
                      "__spacer__",
                      "Standard (clean) line")
  
  method.labels <<- c("Standard",
                      "Standard (clean)",
                      "Adaptive+",
                      "Adaptive+ (asymptotic)",
                      "Adaptive+ (label-cond)",
                      "",
                      "Standard (clean, simple)")
  color.scale <<- cbPalette[c(1,3,4,7,NA,10)]
  shape.scale <<- c(1,2,3,5,NA,NA)
  linetype.scale <<- c(1,1,1,1,0,4)
}


make_figure_601 <- function(exp.num, plot.alpha, plot.data="bigearthnet", plot.guarantee="marginal",
                             plot.contamination="uniform",
                             plot.epsilon=0.1,
                             plot.nu=0,
                             plot.n_train=1000,
                             plot.n_clean,
                             plot.optimistic=TRUE,
                             save_plots=FALSE, reload=FALSE) {
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data==plot.data, n_train==plot.n_train, Guarantee==plot.guarantee,
           Label=="marginal", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, n_clean %in% plot.n_clean,
           n_cal<30000)
  
  df.clean.values <- df %>%
    filter(Method=="Standard using clean") %>%
    group_by(Key) %>%
    summarise(mean_values=mean(Mean))
  df.clean.coverage <- as.numeric(df.clean.values[1,2])
  df.clean.size <- as.numeric(df.clean.values[2,2])
  df.clean <- tibble(Key=c("Coverage","Size"), Mean=c(df.clean.coverage,df.clean.size))
  
  df.clean.legend <- df %>%
    group_by(Key) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal)) %>%
    left_join(df.clean, by = "Key") %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "Standard (clean) line", SE = 0)
  
  df.spacer.legend <- df %>%
    group_by(Key) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal)) %>%
    left_join(df.clean, by = "Key") %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "__spacer__", SE = 0)
  
  # Aggiunge la riga fittizia al df principale (senza Standard using clean)
  df.plot <- df %>%
    filter(Method != "Standard using clean") %>%
    bind_rows(df.clean.legend) %>%
    bind_rows(df.spacer.legend)
  
  df.nominal <- tibble(Key="Coverage", Mean=1-plot.alpha)
  df.range <- tibble(Key=c("Coverage","Coverage"), Mean=c(0.89,0.92), n_cal=1000, Method="Standard")
  
  pp <- df.plot %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    ggplot(aes(x=n_cal, y=Mean, color=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.1) +
    facet_wrap(~Key, scales = "free_y", labeller = as_labeller(c("Coverage" = "Coverage", "Size" = "Size"))) +
    #facet_grid(Key~N_CLEAN, scales="free") +
    geom_hline(data=df.nominal, aes(yintercept=Mean), linetype="dashed") +
    geom_point(data=df.range, aes(x=n_cal, y=Mean), alpha=0) +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    scale_x_continuous(trans='log10') +
    xlab("Number of noisy calibration samples") +
    ylab("") +
    guides(                                                                # <--
      color    = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1))),
      shape    = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1))),
      linetype = guide_legend(override.aes = list(alpha = c(1,1,1,1,0,1)))
    ) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.margin = margin(5, 1, 1, -10))
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_nt%d_ncl%d_eps%s_nu%s_%s_optimistic%s.pdf",
                         exp.num, plot.data, plot.n_train, plot.n_clean, plot.epsilon, plot.nu, plot.contamination, plot.optimistic)
    ggsave(file=plot.file, height=2.5, width=9, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

exp.num <- 601
plot.data <- "bigearthnet"
plot.alpha <- 0.1
plot.epsilon <- 0.016
plot.contamination <- "real"
plot.n_train <- 5000
plot.n_clean <- c(500)
plot.guarantee="marginal"

## Not shown in paper
#make_figure_601(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.guarantee="marginal",
#                 plot.contamination=plot.contamination,
#                 plot.epsilon=plot.epsilon,
#                 plot.n_train=plot.n_train, plot.n_clean=plot.n_clean,
#                 save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)

#### Figure with zoom for the paper -----------------------------------------
library(cowplot)
library(patchwork)
library(ggforce)


get_legend2 <- function(pl) {
  g <- ggplotGrob(pl)
  idx <- grep("guide-box", g$layout$name)
  legends <- g$grobs[idx]
  legends <- legends[!vapply(legends, function(x) inherits(x, "zeroGrob"), logical(1))]
  if (length(legends) == 0) stop("No legend found - check that legend.position != 'none' on the plot you extract from.")
  legends[[1]]
}

make_figure_601b <- function(exp.num, plot.alpha, plot.data="bigearthnet", plot.guarantee="marginal",
                              plot.contamination="uniform",
                              plot.epsilon=0.1,
                              plot.nu=0,
                              plot.n_train=1000,
                              plot.n_clean,
                              plot.optimistic=TRUE,
                              zoom.ylim=c(1.30,1.38),
                              save_plots=FALSE, reload=FALSE) {
  
  if(reload) {
    summary <- load_data(exp.num)
  }
  
  init_settings(plot.optimistic = plot.optimistic)
  
  df <- summary %>%
    filter(data==plot.data, n_train==plot.n_train, Guarantee==plot.guarantee,
           Label=="marginal", Alpha==plot.alpha,
           Method %in% method.values,
           contamination==plot.contamination,
           epsilon==plot.epsilon, n_clean %in% plot.n_clean)
  
  df.clean.values <- df %>%
    filter(Method=="Standard using clean") %>%
    group_by(Key) %>%
    summarise(mean_values=mean(Mean))
  df.clean.coverage <- as.numeric(df.clean.values[1,2])
  df.clean.size     <- as.numeric(df.clean.values[2,2])
  df.clean <- tibble(Key=c("Coverage","Size"), Mean=c(df.clean.coverage,df.clean.size))
  
  df.clean.legend <- df %>%
    group_by(Key) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal)) %>%
    left_join(df.clean, by = "Key") %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "Standard (clean) line", SE = 0)
  
  df.spacer.legend <- df %>%
    group_by(Key) %>%
    summarise(n_cal_min = min(n_cal), n_cal_max = max(n_cal)) %>%
    left_join(df.clean, by = "Key") %>%
    tidyr::pivot_longer(c(n_cal_min, n_cal_max), values_to = "n_cal") %>%
    mutate(Method = "__spacer__", SE = 0)
  
  df.plot <- df %>%
    filter(Method != "Standard using clean") %>%
    bind_rows(df.clean.legend) %>%
    bind_rows(df.spacer.legend) %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels))
  
  df.coverage <- df.plot %>% filter(Key == "Coverage")
  df.size     <- df.plot %>% filter(Key == "Size")
  
  base_theme <-
    theme_bw() +
    theme(
      text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  p.coverage <-
    ggplot(df.coverage, aes(n_cal, Mean, colour=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.1) +
    geom_hline(yintercept=1-plot.alpha, linetype="dashed") +
    scale_x_log10() +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    labs(x="Number of noisy calibration samples", y="Coverage") +
    base_theme
  
  p.size <-
    ggplot(df.size, aes(n_cal, Mean, colour=Method, shape=Method, linetype=Method)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.1) +
    scale_x_log10() +
    scale_color_manual(values=color.scale) +
    scale_shape_manual(values=shape.scale) +
    scale_linetype_manual(values=linetype.scale) +
    labs(x="Number of noisy calibration samples", y="Size") +
    base_theme

  p.size.zoom <- p.size +
    facet_zoom(ylim = zoom.ylim, zoom.size = 2, horizontal = TRUE, show.area = TRUE)
  
  # Build the legend from a throwaway copy of p.size with the legend turned on
  legend_src <- p.size + theme(legend.position = "right", legend.title = element_text(size=12))
  legend <- get_legend2(legend_src)

  figure <-
    (p.coverage | wrap_elements(full = p.size.zoom) | wrap_elements(legend)) +
    plot_layout(widths = c(0.8, 1.05, 0.5))
  
  g <- figure
  
  if(save_plots) {
    plot.file <- sprintf("figures/exp%d_%s_nt%d_ncl%d_eps%s_nu%s_%s_optimistic%s_zoom.pdf",
                         exp.num, plot.data, plot.n_train, plot.n_clean, plot.epsilon, plot.nu, plot.contamination, plot.optimistic)
    ggsave(plot.file, plot=g, height=2.5, width=9, units="in")
    return(NULL)
  } else {
    return(g)
  }
}

exp.num <- 601
plot.data <- "bigearthnet"
plot.alpha <- 0.1
plot.epsilon <- 0.016
plot.contamination <- "real"
plot.n_train <- 5000
plot.n_clean <- c(500)
plot.guarantee <- "marginal"

# Figure 6
make_figure_601b(exp.num=exp.num, plot.alpha=plot.alpha, plot.data=plot.data, plot.guarantee=plot.guarantee,
                  plot.contamination=plot.contamination,
                  plot.epsilon=plot.epsilon,
                  plot.n_train=plot.n_train, plot.n_clean=plot.n_clean,
                  zoom.ylim=c(1.30,1.4),
                  save_plots=TRUE, plot.optimistic=TRUE, reload=TRUE)
dev.off()






# IMAGES APPENDIX B5: Simplified Methods for Special Contamination Models ------

method.values = c("RR", "CVX")
method.labels = c("RR", "Opt")
cbPalette <- c("grey50", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
color.scale <- cbPalette[c(3,7)]
shape.scale <- c(2,0)
linetype.scale <- c(1,1)

load_data <- function() {
  idir <- "results/simplified_methods"
  ifile.list <- list.files(idir, recursive=FALSE)
  ifile.list <- ifile.list[!grepl("comparison", ifile.list)]
  
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  
  summary <- results %>%
    pivot_longer("n_cal", names_to = "Key", values_to = "Value")
  
  return(results)
}

### Block randomized response model  -------------------------------------------
# Plot of the finite sample correction for the optimized betas and the rr betas

#' Figure A1: Finite sample correction as a function of the calibration set size,
#'         for different number of classes
make_figure_A1 <- function(save_plots=FALSE, plot.epsilon=0.1, plot.model_name="B",
                           method.values, method.labels, label.values, label.labels) {

  summary <- load_data()
  
  df <- summary %>%
    filter(plot=="Klab", epsilon==plot.epsilon, model_name==plot.model_name, K %in% label.values) %>%
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
    plot.file <- sprintf("figures/delta_marg_const_%s_Klab_eps%f.pdf", plot.model_name, plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(4, 8, 16)
label.labels = c("K=4", "K=8", "K=16")
plot.model_name = "B"
plot.epsilon = 0.1

# Figure A1
make_figure_A1(plot.epsilon=plot.epsilon, plot.model_name=plot.model_name, method.values=method.values, method.labels=method.labels,
               label.values=label.values, label.labels=label.labels, save_plots=TRUE)

#' Figure A2: Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size
make_figure_A2 <- function(save_plots=FALSE, plot.epsilon=0.1, plot.model_name="B") {

  summary <- load_data()
  
  df <- summary %>%
    filter(plot=="nlab", epsilon==plot.epsilon, model_name==plot.model_name, n_cal %in% label.values) %>%
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
    plot.file <- sprintf("figures/delta_marg_const_%s_nlab_eps%f.pdf", plot.model_name, plot.epsilon)
    ggsave(file=plot.file, height=3, width=7, units="in")
    return(NULL)
  } else{
    return(pp)
  }
}

label.values = c(100, 500, 1000)
label.labels = c("n=100", "n=500", "n=1000")
plot.epsilon = 0.1
plot.model_name = "B"

# Figure A2
make_figure_A2(plot.epsilon=plot.epsilon, plot.model_name=plot.model_name, save_plots=TRUE)

### Two-level randomized response model ----------------------------------------
### Plot of the finite sample correction for the optimized betas and the rr betas

#' Figure A3: Finite sample correction as a function of the calibration set size,
#'         for different number of classes and for different combinations of epsilon and nu
#'

make_figure_A3 <- function(save_plots=FALSE, plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8)) {
  
  summary <- load_data()
  
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

# Figure A3
make_figure_A3(plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8), save_plots=TRUE)



#' Figure A4 : Finite sample correction as a function of the number of classes, for different
#'         values of the calibration set size and for different combinations of epsilon and nu
#'         
make_figure_A4 <- function(save_plots=FALSE, plot.epsilon=0.1, plot.nu.vals=c(0.2,0.8)) {
  
  summary <- load_data()
  
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

make_figure_A4(plot.epsilon=0.1, plot.nu.vals=c(0.2, 0.8), save_plots=TRUE)

#' Figure A5: Finite sample correction as a function of nu, for different
#'         values of the calibration set size and for different number of classes
#'         
make_figure_A5 <- function(save_plots=FALSE, plot.epsilon=0.1, plot.K.vals=c(4,16), model_name) {

  summary <- load_data()
  
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
plot.K.vals = c(4,16)
model_name = "BRR"

# Figure A5
plot.epsilon = 0.1
make_figure_A5(plot.epsilon=plot.epsilon, plot.K.vals=plot.K.vals, model_name=model_name,
               save_plots=TRUE)

# Figure A6
plot.epsilon = 0.5
make_figure_A5(plot.epsilon=plot.epsilon, plot.K.vals=plot.K.vals, model_name=model_name,
               save_plots=TRUE)
