# Particle MCMC algorithm to sample parameters of a velocity jump process model 
# Can be used to examine effects of sampling frequency and measurment noise on inference for a velocity jump process model.
# For further details including how to run with arguments from command line, see README
# Author: Jonathan Harrison, University of Oxford, 11/07/2017

#######################################################
require(MASS,warn.conflicts=FALSE)
require(mcmc,warn.conflicts=FALSE)

args = commandArgs(trailingOnly=TRUE) #read in arguments from the command line
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  identifier = 'simplifiedv01'
  hyperparams = matrix(NA,nrow=0,ncol=0)
  hyperparams['lambda'] = 0.1
  hyperparams['sigma'] = 0.05
  hyperparams['speed'] = 1
  hyperparams['dt'] = 4
  hyperparams['finaltime'] =  2^7
  hyperparams['seed'] = 123  
  hyperparams['M'] = 50000 #numerber of MCMC iterations
  hyperparams['N'] = 400 #numer of particles in the filter
  cat("No inputs supplied. Using default value to look for data to read in.")
} else if (length(args)>=1) {
  identifier = args[1]
  hyperparams = matrix(NA,nrow=0,ncol=0)
  hyperparams['lambda'] = 0.1
  hyperparams['sigma'] = as.numeric(args[2])
  hyperparams['speed'] = 1
  hyperparams['dt'] = as.numeric(args[3])
  hyperparams['finaltime'] = as.numeric(args[4])
  hyperparams['seed'] = as.integer(args[5])  
  hyperparams['M'] = 5000 #numerber of MCMC iterations
  hyperparams['N'] = 400 #numer of particles in the filter
}

ptm <- proc.time()
# Check for dependent files
dependency_check =  file.exists("bootstrap_filter_simplified.r")   &&
  file.exists("emission_prob_simplified.r")              &&
  file.exists("observe_vjump.r")              &&
  file.exists("vjump.r")              &&
  file.exists("reorientation_kernel.r")              &&
  file.exists("pij.cpp")              &&
  file.exists("./data/")              &&
  file.exists("./plots/")             &&
  file.exists("./data/cov_est")
if ( !dependency_check ) {
  stop( "Program requires the files and folders: bootstrap_filter.r, cov_est, data, plots. 
        Perhaps check the working directory is set to the right one?")
}

# Read in data and functions
source("bootstrap_filter_simplified.r")
source("emission_prob_simplified.r")
source('observe_vjump.r')
source('vjump.r')
source('reorientation_kernel.r')
Rcpp::sourceCpp('pij.cpp') # load cpp code

 y <- observe_vjump(hyperparams) #generate dataset  

cat(paste('Simulation params are: \n'))
print(hyperparams)
if (hyperparams["sigma"] > 0) {
  is.data.noisy = TRUE
} else {
  is.data.noisy = FALSE
  #fix numerical rounding errors
  y[abs(y)<10^-12] = 0
}

################################

set.seed(100)

#### Initialisation ####
#set sampling frequency
dt = hyperparams['dt']
# Initialise paramters 
params_cur = c( log10(hyperparams['lambda']),                      
                log10(hyperparams['sigma']))
params_prop = params_cur
# Set number of time steps
t_end = length(y)
if (t_end>10^4) stop('length of time series is too long. t_end>10^4') #run time likely to be very long in this case
# Set number MCMC iterations
M = hyperparams['M']
# Set number particles
N = hyperparams['N']
# Set thinning parameter
nthin = 2
cat(paste('M:',M,'N:',N,'\n',sep=' '))
# Run bootstrap filter to obtain current state
current_state = bootstrap_filter_simplified( t_end, N, 10^(params_prop), y, dt, is.data.noisy)
prop_state = current_state
loglik_cur = current_state$loglik
loglik_prop = loglik_cur
# Set proposal sds for random walk
proposal_sd = as.matrix( read.table( "./data/cov_est", header=TRUE, row.names=1 ) )
# State storage
state_container = matrix(NA, nrow=M, ncol=3+t_end)
colnames(state_container) = c( "loglik", "lambda", "sigma", paste( rep("x", t_end), 1:t_end, sep="" ) )
acc_prob = 0
param_conditions = 0
av_acc = 0
##################################

#### Recursion ####
for ( iter1 in 1:M ) {
  
  for ( iter2 in 1:nthin ) {
    
    # Propose a new state and estimate its marginal likelihood
    params_prop = mvrnorm( 1, params_cur, proposal_sd )
    param_conditions =  (params_prop["lambda"]   > -2) && 
      (params_prop["lambda"] < 1) &&
      (params_prop["sigma"] >=  -Inf) && 
      (params_prop["sigma"] < 1)
    if ( param_conditions == FALSE ) {
      # If the parameter conditions are not met then reject the sample
      acc_prob = 0
    } else {
      # estimate marginal likelihood of proposed state using a particle filter
      prop_state = bootstrap_filter_simplified( t_end, N, 10^(params_prop), y, dt, is.data.noisy)
      loglik_prop = prop_state$loglik
      acc_prob = exp( loglik_prop - loglik_cur )
    }
    # Decide whether to accept state
    if ( runif(1) < acc_prob ) {
      loglik_cur = loglik_prop
      params_cur = params_prop
      current_state = prop_state
      # Calculate average acceptance
      av_acc = av_acc + 1/(M*nthin)
    } 
  }
  if (iter1%%100==0){  #only print every 100 iterations
    cat(paste(iter1," ",sep=""))
  }
  state_container[iter1,] = c( loglik_cur, params_cur, current_state$path )
}
#################################

av_out = colMeans( state_container )
ess_est <- coda::effectiveSize( state_container[,2] )
print(ess_est)
# Output whole sample
write.table( state_container, file=paste('./data/mcmc_',identifier,'.out',sep=''), sep="\t" )
# Ouput mean parameters and hidden trajectory
write.table( av_out, file=paste('./data/mcmc_',identifier,'.mean',sep=''), sep="\t" )
# Output acceptance probability
write.table( av_acc, file=paste('./data/mcmc_',identifier,'.acc',sep=''),sep="\t" )

#### Plot and Store Results ####
# Output traceplots and acf plots for each parameter
pdf(paste('./plots/lambda_ts_',identifier,'.pdf',sep=''))
plot( ts(state_container[,2]), ylab="Lambda", xlab="Iteration", col="maroon" )
dev.off()
pdf(paste('./plots/lambda_acf_',identifier,'.pdf',sep=''))
plot( acf( state_container[,2]), main="" )
dev.off()
if (max(state_container[,3])>-Inf){
  #then sigma was not originally 0 and we have a trace plot for sigma
  pdf(paste('./plots/sigma_ts_',identifier,'.pdf',sep=''))
  plot( ts(state_container[,3]), ylab="Sigma", xlab="Iteration", col="maroon" )
  dev.off()
  pdf(paste('./plots/sigma_acf_',identifier,'.pdf',sep=''))
  plot( acf( state_container[,3]), main="" )
  dev.off() # can uncomment when also running chain for sigma
}
# Output covariance, which could be used for tuning
cov_est <- c(initseq(state_container[,2])$var.con)
write.table(cov_est, file=paste('./data/cov_est_',identifier,'.mean',sep=''),sep="\t" )
lambda.95pc.posterior = quantile(state_container[,2],probs=c(0.025,0.5,0.975))
print(lambda.95pc.posterior)

print(proc.time() - ptm)
