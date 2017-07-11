bootstrap_filter_simplified = function( t_end, N, params, y, dt, is.data.noisy)
{
  # Bootstrap Filter
  # t_end: end time
  # N: number of particles
  # params: vector of parameters  *** do not input log10 of the params ***
  # y: vector of observed angle changes
  # dt: time discretization
  # is.data.noisy: logical for whether measurement noise is present, for calculation of emission probabilities
  
  # x_t ~ Bernoulli(lambda*dt) 
  # y_t | x_t given by emission probabilities derived in manuscript
  
  
  #### Initialise ####
  mean.boot = rep( NA, t_end )
  # path and index containers
  path_mat = matrix( NA, ncol=N, nrow=t_end )
  index_mat = matrix( NA, ncol=N, nrow=t_end-1 )
  # Sample particles for initial state x1
  #particles gives the current hidden state for each of the particles in the filter
  particles = (runif(N) < params["lambda"]*dt)
  #assume start off, no transition in first step (no way to detect this), some probability of trabsition after this
  log.w = emission_prob_simplified(y[1],0,rep(0,N),particles,params,is.data.noisy)
  # Estimate log of likelihood estimate, var loglik and normalise weights (log rather than actual for numerical stability)
  b = max(log.w)
  loglik = b + log( mean(exp(log.w - b)) )
  varloglik = log( var(exp(log.w))/N )
  w = exp(log.w-b)/sum(exp(log.w-b)) #w = w/sum(w)
  if (b < -14) {
    warning('sum of weights is 0, particle filter has bombed')
    return( list( "mean"=( exp( mean.boot ) ), "loglik"=-10^99, "path"=path_mat[,1] ) )
  }
  # Store values
  path_mat[1,] = particles
  mean.boot[1] = sum( w * particles )
#######################################  

  #### Recursion ####
  for ( i in 2:t_end ) { 
    # Resample
    index = sample( 1:N, prob = w, size = N, rep = T )
    index_mat[i-1,] = index
    particles = particles[index]
    # Propagate Particles
    #here only depend on j for transition probabilities A_ij
    update_prob = params["lambda"]*dt  
    particles = (runif(N) < update_prob)  #draw bernoulli RVs with p=dt*lambda 
    path_mat[i,] = particles
    # Update weights
    #!!!!!!!!!!!!!!!!!!!! Measurement model
    log.w =  emission_prob_simplified(y[i],y[i-1],path_mat[i-1,],path_mat[i,],params,is.data.noisy)  #emission probability given i,j,k
    # Update log of likelihood estimate and normalise weights
    b = max(log.w) #using log sum exp trick, see http://timvieira.github.io/blog/post/2014/02/11/exp-normalize-trick/ for details 
    loglik = loglik + b + log( mean(exp(log.w-b)) )
    sum.w = sum(exp(log.w-b))
if (b < -14) {
    warning('sum of weights is 0, particle filter has bombed')
    return( list( "mean"=( exp( mean.boot ) ), "loglik"=-10^99, "path"=path_mat[,1] ) )
}
    w = exp(log.w-b)/sum.w #w = w/sum(w)
    # Record mean
    mean.boot[i] = sum( w * particles )
  }
  #### Pick a random path to return ####
  path = rep( NA, t_end )
  index = sample( 1:N, prob = w, size = 1 )
  path_index = c( index_mat[,index], index )
  for ( i in 1:t_end ) {
    path[i] = path_mat[i,path_index[i]]
  }
#######################################
  return( list( "mean"=( exp( mean.boot ) ), "loglik"=loglik, "path"=path ) )
}
