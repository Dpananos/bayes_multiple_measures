library(rstan)

# y<- c(PP, PN, NP, NN)
dl <- list(y = c(164, 415, 268, 3310), 
          N = 4157,
          p_mu = 0.5, 
          p_kappa = 2)

model <- stan_model(file = 'stan/combined_tests.stan')

fit <- sampling(model, data = dl, chains = 12)

sink('fits/combined.txt')
fit
sink()
