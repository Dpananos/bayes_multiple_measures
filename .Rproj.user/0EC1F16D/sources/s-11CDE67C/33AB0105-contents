library(rstan)
library(parallel)
library(tidyverse)
library(abind)
library(ggridges)
source('functions.R')

NUMDRAWS = 20000

#means and standard errors from sampling distribution of means and variances.
# See table 2
#Survey
S2.mean = 0.553
S2.sd = (0.686 - S2.mean)/1.96

C2.mean = 0.937
C2.sd = (0.974 - C2.mean)/1.96

#Admin
S1.mean = 0.629
S1.sd = (0.668 - S1.mean)/1.96

C1.mean = 0.938
C1.sd = (0.947 - C1.mean)/1.96

#Double check moment matching works
##Survey
ggplot()+
  stat_function(aes(x = seq(0,1,0.01)),xlim = c(0,1), n = 1001, fun = function(x) dnorm(x, S2.mean, S2.sd))+
  stat_function(aes(x = seq(0,1,0.01)),xlim = c(0,1), n = 1001, fun = function(x) dbeta(x, return_a(S2.mean, S2.sd), return_b(S2.mean, S2.sd)), color = 'red' )

ggplot()+
  stat_function(aes(x = seq(0,1,0.01)),xlim = c(0,1), n = 1001, fun = function(x) dnorm(x, C2.mean, C2.sd))+
  stat_function(aes(x = seq(0,1,0.01)),xlim = c(0,1), n = 1001, fun = function(x) dbeta(x, return_a(C2.mean, C2.sd), return_b(C2.mean, C2.sd)), color = 'red' )


##Admin
ggplot()+
  stat_function(aes(x = seq(0,1,0.01)),xlim = c(0,1), n = 1001, fun = function(x) dnorm(x, S1.mean, S1.sd))+
  stat_function(aes(x = seq(0,1,0.01)),xlim = c(0,1), n = 1001, fun = function(x) dbeta(x, return_a(S1.mean, S1.sd), return_b(S1.mean, S1.sd)), color = 'red' )

ggplot()+
  stat_function(aes(x = seq(0,1,0.01)),xlim = c(0,1), n = 1001, fun = function(x) dnorm(x, C1.mean, C1.sd))+
  stat_function(aes(x = seq(0,1,0.01)),xlim = c(0,1), n = 1001, fun = function(x) dbeta(x, return_a(C1.mean, C1.sd), return_b(C1.mean, C1.sd)), color = 'red' )


#Do the single simulations.
#----Survey----
survey.prior.params = list(
  'Pi.prior.a' = 1, 'Pi.prior.b' = 1,
  'S.prior.a' = return_a(S2.mean, S2.sd), 'S.prior.b' = return_b(S2.mean, S2.sd),
  'C.prior.a' = return_a(C2.mean, C2.sd), 'C.prior.b' = return_b(C2.mean, C2.sd))

#Detect number of cores on the machine.
numCores<- detectCores()
#Create a cluster locally.
cl <- makeCluster(numCores)
clusterExport(cl,list('survey.prior.params','single.diagnostic.simulation','NUMDRAWS'))

survey.results = parLapply(cl,
                    1:numCores,
                    function(x) single.diagnostic.simulation(a = 579,
                                                             b =3578,
                                                             prior.params = survey.prior.params,
                                                             chainID = x,
                                                             N.draws = NUMDRAWS))

survey.draws = map_dfr(survey.results, bind_rows) %>% 
              group_by(chain) %>% 
              slice((NUMDRAWS/2):NUMDRAWS) %>% 
              ungroup
#Kill the cluster
stopCluster(cl)

#----Admin----
admin.prior.params = list(
  'Pi.prior.a' = 1, 'Pi.prior.b' = 1,
  'S.prior.a' = return_a(S1.mean, S1.sd), 'S.prior.b' = return_b(S1.mean, S1.sd),
  'C.prior.a' = return_a(C1.mean, C1.sd), 'C.prior.b' = return_b(C1.mean, C1.sd)
)


#Detect number of cores on the machine.
numCores<- detectCores()
#Create a cluster locally.
cl <- makeCluster(numCores)
clusterExport(cl,list('admin.prior.params','single.diagnostic.simulation','NUMDRAWS'))

admin.results = parLapply(cl,
                           1:numCores,
                           function(x) single.diagnostic.simulation(a = 432,
                                                                    b = 3725,
                                                                    prior.params = admin.prior.params,
                                                                    chainID = x,
                                                                    N.draws = NUMDRAWS) )

admin.draws = map_dfr(admin.results, bind_rows) %>% 
              group_by(chain) %>% 
              slice((NUMDRAWS/2):NUMDRAWS) %>% 
              ungroup
#Kill the cluster
stopCluster(cl)


#----Both----

prior.params = list(
  #----Prevalence priors----
  'Pi.prior.a' = 1,
  'Pi.prior.b' = 1,
  'S2.prior.a' = return_a(S2.mean, S2.sd), 'S2.prior.b' = return_b(S2.mean, S2.sd),
  'C2.prior.a' = return_a(C2.mean, C2.sd), 'C2.prior.b' = return_b(C2.mean, C2.sd),
  'S1.prior.a' = return_a(S1.mean, S1.sd), 'S1.prior.b' = return_b(S1.mean, S1.sd),
  'C1.prior.a' = return_a(C1.mean, C1.sd), 'C1.prior.b' = return_b(C1.mean, C1.sd)
)

#Detect number of cores on the machine.
numCores<- detectCores()
#Create a cluster locally.
cl <- makeCluster(numCores)
clusterExport(cl,list('prior.params','two.diagnostics.simulation','NUMDRAWS'))

dual.results = parLapply(cl,
                    1:numCores,
                    function(s) two.diagnostics.simulation(
                      u = 164,
                      v = 415,
                      w = 268,
                      x = 3310,
                      prior.params =prior.params, 
                      chainID = s, 
                      N.draws = NUMDRAWS))

#Kill cluster

stopCluster(cl)
dual.draws = map_dfr(dual.results, bind_rows) %>% 
              group_by(chain) %>% 
              slice((NUMDRAWS/2):NUMDRAWS) %>% 
              ungroup

#----Convergence----
check.convergence(admin.draws)
check.convergence(survey.draws)
check.convergence(dual.draws)



ggplot()+
  geom_density(data = dual.draws , aes(Pi), color = 'black', fill = 'black', alpha = 0.25)+
  geom_density(data = survey.draws , aes(Pi), color = 'red', fill = 'red', alpha = 0.25)+
  geom_density(data = admin.draws , aes(Pi),color = 'blue', fill = 'blue', alpha = 0.25)+
  scale_x_continuous(limits = c(0,.3))+
  theme_minimal()+
  labs(x = 'Prevalence', y = 'Density', fill="Measure")+
  theme(aspect.ratio = 1/1.61)




ggplot()+
  stat_summary(data = dual.draws, aes('Dual',Pi), fun.data = mean_qi)+
  stat_summary(data = survey.draws, aes('Survey',Pi), fun.data = mean_qi)+
  stat_summary(data = admin.draws, aes('Admin',Pi), fun.data = mean_qi)+
  scale_y_continuous(labels = scales::percent)+
  labs(x = 'Instrument', y = 'Estimated Prevelance')+
  coord_flip()+
  theme(aspect.ratio = 1/1.61)


pis = tibble(instrument = 'dual', pi = dual.draws$Pi) %>% 
      bind_rows(tibble(instrument = 'survey', pi = survey.draws$Pi)) %>% 
      bind_rows(tibble(instrument = 'admin', pi = admin.draws$Pi)) 


pis %>% 
  ggplot(aes(y = instrument, x = pi))+
  geom_density_ridges2(rel_min_height = 0.01, scale = 0.75)+
  xlim(0,0.25)+
  labs(y = 'Instrument', x = 'Estimated Prevelance')+
  theme_minimal()


fig_data = pis %>% 
  mutate(instrument = stringr::str_to_title(instrument),
         instrument = if_else(instrument=='Dual','Combined',instrument),
         instrument = factor(instrument, levels = c('Admin','Survey','Combined'), ordered = T))


theme_set(theme_bw())
fig_data %>% 
  ggplot(aes(x = pi, fill = instrument))+
  geom_histogram(aes(y=..density..),color = 'black', binwidth = 0.01)+
  facet_wrap(~instrument, nrow = 3)+
  scale_x_continuous(labels = scales::percent)+
  scale_fill_brewer(palette = 'Set1')+
  labs(x = 'Posterior Population Prevalence',
       y = '')+
  guides(fill = F)

ggsave('fig1_1.png', dpi = 400, height = 5, width = 5)




fig_data %>% 
  ggplot(aes(x = instrument, y = pi, color = instrument))+
  stat_summary(fun.data = tidybayes::mean_qi, )+
  scale_y_continuous(labels = scales::percent)+
  scale_color_brewer(palette = 'Set1')+
  labs(x = 'Instrument',
       y = 'Posterior Population Prevalence'
       )+
  guides(color = F)+
  coord_flip()+
  theme(aspect.ratio = 1/1.61)

ggsave('fig1_2.png', dpi = 400, height = 5, width = 5)

