library(rstan)
library(tidyverse)
library(tidybayes)
library(cowplot)
library(posterior)

theme_set(theme_classic())

# y<- c(PP, PN, NP, NN)
dl <- list(y = c(164, 415, 268, 3310), 
          N = 4157,
          p_mu = 0.5, 
          p_kappa = 2,
          mu_se_1 = 0.692,
          sd_se_1 = 0.02,
          mu_se_2 = 0.553, 
          sd_se_2 = 0.07,
          mu_sp_1 = 0.938,
          sd_sp_1 = 0.005,
          mu_sp_2 = 0.937,
          sd_sp_2 = 0.02)

model <- stan_model('stan/combined_tests.stan')

fit <- sampling(model, data = dl, chains = 12)

sink('fits/combined.txt')
print(fit, digits=3, probs = c(0.025, 0.5, 0.975))
sink()

cells = tibble(i = factor(1:4), 
               cells = c('P/P','N/P','P/N','N/N'),
               vals = c(164, 415, 268, 3310))

fit %>% 
  gather_draws(y_ppc[i]) %>% 
  mutate(i = factor(i)) %>% 
  left_join(cells) %>% 
  ggplot()+
  stat_histinterval(aes(.value, cells),size = 1)+
  geom_point(aes(vals, cells), color = 'red', size = 0.5)+
  theme(aspect.ratio = 1/1.61)+
  labs(x='Estimated cell count',
       y='Cell')
  
ggsave('figures/predictive_check.png', width = 3*1.61, height = 3)


