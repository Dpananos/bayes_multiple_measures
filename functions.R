#Need these to convert prior studies into Beta parameters
return_a <- function(mu, sigma){
  x = ((1-mu)/sigma^2 - 1/mu)*mu^2
  return(x)
}

return_b <- function(mu,sigma){
  x = return_a(mu,sigma)*(1/mu - 1)
  return(x)
}



two.diagnostics.simulation<-function(u,v,w,x,prior.params, chainID, N.draws = 25000, seed = NULL){
  if(any(c(u,v,w,x)<0)) stop('Contignency table frequencies must be non-negative')
  if(N.draws<0) stop("Are you dumb?!")
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  #Because there are a lot of parameters to pass
  #I find it easier to put them in a named list
  #and then pass the named list.
  #If there are any errors, this would be the place to look for them
  env = environment()
  list2env(prior.params, envir = env)
  
  
  #Initialize vector for samples to live in.
  #Then draw a seed sample from the prior.
  Pi.samples<-rep(0,N.draws)
  Pi<-rbeta(1,Pi.prior.a, Pi.prior.b)
  
  S1.samples<-rep(0, N.draws)
  S1<-rbeta(1,S1.prior.a, S1.prior.b)
  
  S2.samples<-rep(0, N.draws)
  S2<-rbeta(1,S2.prior.a, S2.prior.b)
  
  C1.samples<-rep(0,N.draws)
  C1<-rbeta(1,C1.prior.a, C1.prior.b)
  
  C2.samples<-rep(0,N.draws)
  C2<-rbeta(1,C2.prior.a, C2.prior.b)
  
  #Initalize a place for the pseudo data to live
  Y1.samples<-rep(0,N.draws)
  Y2.samples<-rep(0, N.draws)
  Y3.samples<-rep(0, N.draws)
  Y4.samples<-rep(0, N.draws)
  
  for( i in 1:N.draws){

    #If we knew the prevalence, sensitivity, and specificity
    #Then the data would look like
    Y1<- rbinom(1,u,Pi*S1*S2/( Pi*S1*S2 + (1-Pi)*(1-C1)*(1-C2)) )
    Y1.samples[i]<- Y1
    
    Y2<- rbinom(1,v ,Pi*S1*(1-S2)/(Pi*S1*(1-S2) + (1-Pi)*(1-C1)*C2))
    Y2.samples[i]<- Y2
    
    Y3<- rbinom(1,w, Pi*(1-S1)*S2/(Pi*(1-S1)*S2 + (1-Pi)*C1*(1-C2)) )
    Y3.samples[i]<-Y3
    
    Y4<- rbinom(1,x, Pi*(1-S1)*(1-S2)/(Pi*(1-S1)*(1-S2) + (1-Pi)*C1*C2 ) )
    Y4.samples[i]<-Y4
    
    
    #If we knew the data,
    #The prevalence, sensitivity, and specificity would look like
    
    #Draw from posterior.  Will be used in next round of sampling.
    Pi<- rbeta(1,Y1+Y2 + Y3 + Y4 + Pi.prior.a, u + v + w + x - Y1 - Y2 - Y3 - Y4 + Pi.prior.b)
    Pi.samples[i]<-Pi
    
    S1<- rbeta(1,Y1 + Y2 + S1.prior.a, Y3+Y4 + S1.prior.b)
    S1.samples[i]<-S1
    
    C1<- rbeta(1,w + x - Y3 - Y4 + C1.prior.a, u+v - Y1 - Y2 + C1.prior.b)
    C1.samples[i]<-C1

    S2<- rbeta(1,Y1+Y3+S2.prior.a, Y2+Y4 + S2.prior.b)
    S2.samples[i] <- S2
    
    C2<- rbeta(1,v+x-(Y2+Y4) + C2.prior.a, u+w - Y1 - Y3 + C2.prior.b)
    C2.samples[i]<- C2
  }
  
  #Return the results as a dataframe
  #TODO:  Might be easier to consider reshaPing so
  #output can be passed directly to stan::monitor.
  results=data.frame(
    'Pi' = Pi.samples,
    'S1' = S1.samples,
    'S2' = S2.samples,
    'C1' = C1.samples,
    'C2' = C2.samples,
    'Y1'= Y1.samples,
    'Y2' = Y2.samples,
    'Y3' = Y3.samples,
    'Y4' = Y4.samples,
    'chain' = chainID,
    'iteration' = 1:N.draws
  ) 
  
  return(results)
  
}



single.diagnostic.simulation<-function(a,b, prior.params, chainID, N.draws = 25000, seed = NULL){
  if(any(c(a,b)<0)) stop('Frequencies must be non-negative')
  if(N.draws<0) stop("Are you dumb?!")
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  #Because there are a lot of parameters to pass
  #I find it easier to put them in a named list
  #and then pass the named list.
  #If there are any errors, this would be the place to look for them
  env = environment()
  list2env(prior.params, envir = env)
  
  
  #Initialize vector for samples to live in.
  #Then draw a seed sample from the prior.
  Pi.samples<-rep(0,N.draws)
  Pi<-rbeta(1,Pi.prior.a, Pi.prior.b)
  
  S.samples<-rep(0, N.draws)
  S<-rbeta(1,S.prior.a, S.prior.b)
  
  C.samples<-rep(0,N.draws)
  C<-rbeta(1,C.prior.a, C.prior.b)
  
  #Initalize a place for the pseudo data to live
  Y1.samples<-rep(0,N.draws)
  Y2.samples<-rep(0, N.draws)
  
  
  for( i in 1:N.draws){
    
    #If we knew the prevalence, sensitivity, and specificity
    #Then the data would look like
    Y1<- rbinom(1,a, Pi*S/(Pi*S + (1-Pi)*(1-C)) )
    Y1.samples[i]<- Y1
    
    Y2<- rbinom(1,b ,Pi*(1-S)/(Pi*(1-S) + (1-Pi)*C) )
    Y2.samples[i]<- Y2
    
    
    #If we knew the data,
    #The prevalence, sensitivity, and specificity would look like
    
    #Draw from posterior.  Will be used in next round of sampling.
    Pi<- rbeta(1,Y1 + Y2 + Pi.prior.a, a+b - Y1 - Y2 + Pi.prior.b)
    Pi.samples[i]<-Pi
    
    S<- rbeta(1,Y1 + S.prior.a, Y2 + S.prior.b)
    S.samples[i]<-S
    
    C<- rbeta(1,b-Y2 + C.prior.a, a-Y1+ C.prior.b)
    C.samples[i]<-C
  }
  
  #Return the results as a dataframe
  #TODO:  Might be easier to consider reshaPing so
  #output can be passed directly to stan::monitor.
  results=data.frame(
    'Pi' = Pi.samples,
    'S' = S.samples,
    'C' = C.samples,
    'chain' = chainID,
    'iteration' = 1:N.draws
  ) 
  
  return(results)
}



check.convergence<-function(draws){
  
  o.draws = draws%>%
    as_tibble() %>% 
    gather(key, var, -iteration, -chain) %>% 
    spread(chain, var)
  
  
  #Monitoring
  vars = o.draws$key
  o.draws = o.draws %>% select(-key, -iteration)
  for.monitoring = split(o.draws, vars )
  admin.monitor.me = abind(for.monitoring, along = 3)
  monitor(admin.monitor.me, 
          warmup = 0, 
          digits_summary = 5, 
          probs = c(0.025,0.5, 0.975),
          se =F)
  return(0)
  
}


mean_qi = function(x){
  
  d = tibble(y = mean(x),
             ymin = quantile(x, 0.025),
             ymax = quantile(x, 0.975))
  
  return(d)
}
