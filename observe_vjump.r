observe_vjump <- function(hyperparams){
#observe a simulated velocity jump process at discrete times
#discrete observations can be misleading, as the apparent direction of travel will not in general be the same as the true direction of travel
#measurement noise may also be present
# depends on vjump.r, reorientation_kernel.r
###############################
require(dplyr,warn.conflicts=FALSE)

if (length(hyperparams)<=0) {  
  hyperparams=matrix(NA,ncol=0,nrow=0)
  hyperparams['lambda'] = 0.1
  hyperparams['sigma'] = 0
  hyperparams['speed'] = 1
  hyperparams['dt'] = 4
  hyperparams['finaltime'] = 2^7
  hyperparams['seed'] = 999  
  cat("No inputs supplied. Using default values.")
}
set.seed(hyperparams['seed'])
output = vjump(hyperparams) #run velocity jump process
theta=output$theta #true angle changes
time=output$time #transition times

pos <- get_true_positions_from_angles(time,theta,hyperparams)
observed.angle.change <- process_raw_path(time,theta,pos[,1],pos[,2],hyperparams)

return(observed.angle.change)
}

#################################################################################
get_true_positions_from_angles <- function(time,theta,hyperparams){
  diff.time <- time[-1] - time[-length(time)]
  pos = matrix(NA,nrow=length(time),ncol=2)
  pos[1,] = c(0,0)
  current_direction = vector('double',length(time))
  current_direction[1] = 0
  for (k in 2:length(time)){
    pos[k,] = pos[k-1,] + hyperparams['speed']*diff.time[k-1]*c(cos(current_direction[k-1]),sin(current_direction[k-1]))
    current_direction[k] = modulo.2pi(current_direction[k-1]+theta[k])
  }
  return(pos)
}

modulo.2pi <- function(x){
  y = x-2*pi*floor(x/(2*pi)+0.5);
  return(y)  
}

process_raw_path <- function(time,theta,posx,posy,hyperparams){
  final.time <- hyperparams['finaltime']
  dt <- hyperparams['dt']
  nt <- 1+final.time/dt #number of time points
  observation.times <- seq(from=0,to=final.time,by=dt)
  interpolatedx <- approxfun(time,posx)
  interpolatedy <- approxfun(time,posy)
  observation.pos <- cbind(interpolatedx(observation.times),interpolatedy(observation.times)) 
  pos.change <- observation.pos[2:nt,] - observation.pos[1:(nt-1),]
  angle.change <- atan2(pos.change[,2],pos.change[,1]) #framewise angle change
  observed.angle.change <- angle.change[2:(nt-1)] - angle.change[1:(nt-2)] + rnorm((nt-2),mean=0,sd=hyperparams['sigma'])
  observed.angle.change = modulo.2pi(observed.angle.change)
  if (near(hyperparams['sigma'],0)) {
    observed.angle.change[abs(observed.angle.change)<10^-12] = 0
  }
  return(observed.angle.change)
}



