library(magrittr) # for %$% operator 
library(tidytext)  # for easy handling of text
library(ldatuning) # for mathematical hint of number of topics
library(topicmodels) # nomen est omen
library(tidyverse)

loadLDAparams <- function(out){
  control_list_gibbs <- list(
    burnin = 250, #2500,
    iter = 500, #5000,
    seed = 0:4,
    nstart = 5,
    best = TRUE
  )
  
  control_list_ctm <- list(
    seed = 0:4,
    nstart = 5,
    best = TRUE
  )
  if(out == "lda") {
    return(control_list_gibbs)
  } else {
    return(control_list_ctm)
  }
  
}


plotGamma <- function(tidymodel) {
  tidymodel %>%
    group_by(document) %>% 
    arrange(desc(gamma)) %>% 
    slice(1) %>% 
    ungroup() %>%
    ggplot(aes(x=gamma)) +
    geom_histogram(bins = 20) +
    xlab("maximum gamma per document") + 
    geom_vline(aes(xintercept = 1/20),
               color="darkred")
  
}

tidy_ctm_gamma  <- function(CTM_object){
  CTM_object %>% 
    slot("gamma")  %>% 
    as_data_frame()  %>% 
    mutate (document = row_number()) %>% 
    gather(topic, gamma, -document) %>%
    mutate(topic = strtoi(stringr::str_sub(topic,2)))
}

tidy_ctm_beta  <- function(CTM_object){
  Terms  <- CTM_object %>% 
    slot("terms") 
  CTM_object %>% 
    slot("beta")  %>% 
    as_data_frame() %>%
    setNames(Terms) %>%
    mutate (topic = row_number()) %>% 
    gather(term, beta, -topic) %>%
    mutate(beta = exp(beta))
}

## helper function to sort 
## correlation matrix
reorder_cor <- function(x){
  ord <- corrplot::corrMatOrder(x)
  x[ord,ord]
}

## helper function to extract 
## lower triangle of matrix
get_lower_tri<-function(x){
  x[upper.tri(x)] <- NA
  return(x)
}

getTopicCorrelations <- function(ctm_gamma, ntopics, plot = TRUE) {
  ## get corr matrix (cor between topics)
  out <- ctm_gamma %>%
    pivot_wider(names_from = topic, values_from = gamma) %>%
    select(-document) %>%
    cor() %>%
    reorder_cor() %>%
    get_lower_tri() %>%
    as_data_frame() %>%
    mutate(topic1 = forcats::as_factor(paste(names(.)))) %>%
    gather(topic2, correlation, - topic1) %>%
    mutate(topic2 = factor(topic2, levels=levels(topic1)))
  
  ## plot corr matrix
  if(isTRUE(plot)) {
    p <- ggplot(out, aes(as.numeric(topic1), as.numeric(topic2), fill=correlation)) +
      geom_tile(color="white") +
      scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4", 
                           midpoint = 0, limit = c(-1,1), space = "Lab",
                           na.value = "white",
                           name="Pearson\nCorrelation") +
      scale_x_continuous(
        breaks=1:ntopics, labels = levels(cor_data$topic1), name="")+
      scale_y_continuous(
        breaks=1:ntopics, labels = levels(cor_data$topic2), name="")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      coord_fixed() 
  }
  if (isTRUE(plot)) {
    return(list(out, p))
  } else {
    return(out)
  }
}
