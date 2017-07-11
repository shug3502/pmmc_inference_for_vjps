vjump = function(hyperparams){
#simulate a path from the velocity jump process model
#essentially uses a Gillespie simulation, ie. a Markov jump process, on the velocity
#and position is updated correspondingly

##################################################

finaltime = hyperparams['finaltime']
lambda = hyperparams['lambda']
maxturns = ceiling(10*finaltime*lambda) #path will be of unknown length, preallocate this much memory
pos = matrix(NA,nrow=maxturns,ncol=2) #store position
pos[1,] = c(0,0) #assume wlog x_0 is at the origin
time = matrix(NA,nrow=maxturns,ncol=1) #store time
time[1]=0
theta=matrix(NA,nrow=maxturns,ncol=1) #store angle changes
theta[1]=0
currenttime=0
k=1

while (currenttime<finaltime){
  k=k+1 #counts number of turns
  tau = rexp(1,rate=lambda) #exponential rv
  theta[k] = reorientation_kernel(theta[k-1])
  pos[k,] = pos[k-1,] + hyperparams['speed']*tau*c(cos(theta[k-1]),sin(theta[k-1]))
  currenttime = currenttime + tau
  time[k] = currenttime
}
posx = pos[1:k,1]
posy = pos[1:k,2]
theta=theta[1:k]
time = time[1:k]
return(list(posx=posx,posy=posy,theta=theta,time=time))
}
