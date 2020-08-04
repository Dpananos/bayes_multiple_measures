library(rstan)
library(tidybayes)
library(tidyverse)
library(posterior)
# Survey Data
# Positive cases, total sample, and prior parameters
admin_dl <- list(y_sample = 432,
           n_sample = 3725+432, 
           spec_mu = 0.938, 
           spec_sd = 0.004, 
           sens_mu = 0.629, 
           sens_sd = 0.020,
           p_mu=0.5,
           p_kappa = 2)

survey_dl <- list(y_sample = 579,
           n_sample = 4157, 
           spec_mu = 0.937, 
           spec_sd = 0.067, 
           sens_mu = 0.553, 
           sens_sd = 0.020,
           p_mu = 0.5, 
           p_kappa = 2)


# Compile the model
model <- stan_model('stan/single_test.stan')

# Sample the model
survey_fit <- sampling(model, data = survey_dl, chains=12)
admin_fit <- sampling(model, data = admin_dl, chains = 12)

sink('fits/survey.txt')
print(survey_fit, digits=3, probs = c(0.025, 0.5, 0.975))
sink()

sink('fits/admin.txt')
print(admin_fit, digits=3, probs = c(0.025, 0.5, 0.975))
sink()

admin_p = admin_fit$draws('p') %>% 
  as_draws_df() %>% 
  spread_draws(p) %>% 
  mutate(model='admin')

survey_p = survey_fit$draws('p') %>% 
  as_draws_df() %>% 
  spread_draws(p) %>% 
  mutate(model='survey')


bind_rows(admin_p, survey_p) %>% 
  ggplot(aes(p, model))+
  stat_histinterval(slab_color = "gray45", 
                    outline_bars = TRUE,
                    breaks = 20)
