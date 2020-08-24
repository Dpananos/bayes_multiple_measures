library(rstan)
library(tidyverse)
library(tidybayes)
library(posterior)
library(magrittr)

param_grid <- crossing(y = list(c(164, 415, 268, 3310)), 
           N = 4157,
           p_mu = 0.5, 
           p_kappa = 2,
           mu_se_1 = seq(0.65, 0.75, 0.05),
           sd_se_1 = 0.02,
           mu_se_2 = seq(0.5, 0.6, 0.05), 
           sd_se_2 = 0.07,
           mu_sp_1 = seq(0.88, 0.98, 0.05),
           sd_sp_1 = 0.005,
           mu_sp_2 = seq(0.88, 0.98, 0.025),
           sd_sp_2 = 0.02) 
  

dl <- transpose(param_grid)
  
model <- stan_model(file = 'stan/combined_tests.stan')

fit_model<-function(dl){
  fit <- sampling(model, data = dl, chains=8, pars = c('p'), include = TRUE)
  
  fit_summ<- fit %>% 
             spread_draws(p) %>% 
             mean_qi()
  
  
  rhats<- bayesplot::rhat(fit)
  divergences<- rstan::get_divergent_iterations(fit)
  fit_summ$rhat <- all(rhats<1.01)
  fit_summ$divs <- any(divergences)
  
  fit_summ
}

sa <- map_dfr(dl, fit_model)

write_csv(sa, 'sensitivity/sensitivity_results.csv')


sa %>% 
  bind_cols(param_grid) %>% 
  mutate(`delta[1]` = mu_se_1, `delta[2]` = mu_se_2) %>% 
  ggplot(aes(mu_sp_2, p, fill = factor(mu_sp_1), color = factor(mu_sp_1)))+
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5)+
  geom_line()+
  facet_grid(`delta[2]` ~ `delta[1]`,
             labeller = function(x) label_parsed(label_both(x, sep= '  == ' )))+
  labs(x=expression(gamma[1]),
       color = expression(gamma[2]),
       fill = expression(gamma[2]),
       y = expression(pi))+
  theme_classic()
  
