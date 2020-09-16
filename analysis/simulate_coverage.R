library(rstan)
library(tidybayes)
library(tidyverse)

# For reproducibility
set.seed(0)

# Compile model
model <- stan_model('stan/combined_tests.stan')

# This function simulates data we need in order to analyze the coverage of our method.
# We draw sensitivities, specificities, and prevalence from our priors
# We run our model and determine how many times our 80% credible interval captures the true prevalence.
make_data<-function(){
  #Draw prevalence
  p <- runif(1)
  
  #Draw a sensitivity and specificty from our priors.
  se_1 <- rnorm(1, 0.629, 0.02)
  se_2 <- rnorm(1, 0.553, 0.07)
  
  sp_1 <- rnorm(1, 0.938, 0.005)
  sp_2 <- rnorm(1, 0.937, 0.02)
  
  #Generate data.
  SE1 <- c(se_1, 1-se_1,se_1, 1-se_1)
  SE2 <- c(se_2, se_2, 1-se_2, 1-se_2)
  SP1 <- c(1-sp_1, sp_1, 1-sp_1, sp_1)
  SP2 <- c(1-sp_2, 1-sp_2, sp_2, sp_2)
  
  p_sample <- p * SE1 * SE2 + (1-p) * SP1 * SP2
  
  y <- as.numeric(rmultinom(1, size = 4157, prob = p_sample))
  
  #Data we need to pass to the model
  dl <- list(y = y, 
             N = 4157,
             p_mu = 0.5, 
             p_kappa = 2,
             mu_se_1 = 0.629,
             sd_se_1 = 0.02,
             mu_se_2 = 0.553, 
             sd_se_2 = 0.07,
             mu_sp_1 = 0.938,
             sd_sp_1 = 0.005,
             mu_sp_2 = 0.937,
             sd_sp_2 = 0.02)
  
  
  list(model_data = dl, p_true = p)
}

extract<-function(data){
  
  #Fit the model
  fit<-sampling(model, data$model_data, chains=4)
  
  #Extract prevalence samples
  p<- fit %>% spread_draws(p) %>% pull(p)
  
  #diagnostics
  rhats<- bayesplot::rhat(fit)
  divergences<- rstan::get_divergent_iterations(fit)
  
  r <- mean_qi(p, .width = c(0.8, 0.9, 0.95))
  
  r %>% 
    mutate(
         divs = sum(divergences),
         rhat = all(rhats<1.01),
         p_true = data$p_true)
}

# Simulate data and fit model 1000 times
results <- rerun(1000, make_data()) %>% 
  map_dfr(extract, .id = 'simulation')

#Save data for later.
write_csv(results, 'coverage/results.csv')

#Plot true prevalence against estimated
read_csv('coverage/results.csv') %>% 
  filter(divs<1, rhat==T) %>% 
  ggplot()+
  geom_pointrange(aes(p_true,y, ymin = ymin, ymax = ymax), color = 'gray')+
  geom_abline()+
  facet_wrap(~.width)
  
#Coverage should be approximately 0.8 since we used an 80% credible interval.
read_csv('coverage/results.csv') %>% 
  group_by(.width) %>% 
  mutate(contained = (p_true<ymax)&(ymin<p_true)) %>% 
  summarise(coverage = mean(contained))
