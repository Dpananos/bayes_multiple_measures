library(tidyverse)
library(rstan)
library(tidybayes)
library(cowplot)
library(posterior)
library(patchwork)

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

# To be used in plotting later
cells <- tibble(i = 1:4,
               counts = dl$y,
               Admin = c('Positive','Negative','Positive','Negative'),
               Survey= c('Positive','Positive','Negative','Negative'))

# Compile the model
model <- stan_model('stan/combined_tests.stan')

# Fit the model
fit <- sampling(model, data = dl, chains = 12)

# Write to a little text file for easy reference later
sink('fits/combined.txt')
print(fit, digits=3, probs = c(0.025, 0.5, 0.975))
sink()

# Make the posterior predictive check plot

make_plot<-function(fit,j){
  
  plotme<-fit %>% 
    gather_draws(y_ppc[i]) %>% 
    left_join(cells) %>% 
    filter(i==j) %>% 
    ggplot()+
    stat_histinterval(aes(.value),
                      outline_bars = T,
                      slab_color = 'gray45')+
    geom_vline(aes(xintercept = counts))+
    labs(x='', y='')+
    theme(aspect.ratio = 1,
          plot.margin = unit(c(1,1,1,1)*0.05, "cm"))
  
  if(j==1){
    plotme<- plotme + 
      facet_grid(~Admin, labeller = function(x) label_both(x,sep = ' ')) + 
      ylab('Frequency')
  }
  if(j==2){
    plotme<- plotme + facet_grid(Survey~Admin, labeller = function(x) label_both(x,sep = ' '))
  }
  if(j==3){
    plotme<-plotme + labs(x='Cell Count', y='Frequency')
  }
  if(j==4){
    plotme<- plotme + 
      facet_grid(Survey~., labeller = function(x) label_both(x,sep = ' ')) + 
      xlab('Cell Count')
  }
  
  plotme
}

pp <- make_plot(fit,1)
np <- make_plot(fit, 2)
pn <- make_plot(fit,3)
nn <- make_plot(fit,4)

final <- (pp + np)/(pn + nn)

ggsave('figures/predictive_check.png', height = 5, width = 5)

#---- Make 6x2 Figure For Paper ----

p = rstan::extract(fit)


# Create subplots
tiff("figures/posterior.tiff", units="in", width=8, height=5, res=300)
par(mfrow=c(2,3))
hist(p$se_1, 
     xlab = expression(delta[1]),
     main = 'Posterior Admin Sensitivity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

hist(p$sp_1, 
     xlab = expression(gamma[1]),
     main = 'Posterior Admin Specificity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

hist(p$p, 
     xlab = expression(pi),
     main = 'Posterior Prevalence',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)


hist(p$se_2, 
     xlab = expression(delta[2]),
     main = 'Posterior Survey Sensitivity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

hist(p$sp_2, 
     xlab = expression(gamma[2]),
     main = 'Posterior Survey Specificity',
     ylab = '',
     yaxt = 'n',
     breaks = 20,
     cex.lab = 1.25)

dev.off()
