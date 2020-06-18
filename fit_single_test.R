library(rstan)

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
model <- stan_model(file = 'stan/single_test.stan')

# Sample the model
survey_fit <- sampling(model, data = survey_dl, chains = 12)
admin_fit <- sampling(model, data = admin_dl, chains = 12)

sink('fits/survey.txt')
survey_fit
sink()

sink('fits/admin.txt')
admin_fit
sink()
