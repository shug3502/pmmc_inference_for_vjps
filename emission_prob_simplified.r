emission_prob_simplified <- function(angle_change, previous_angle_change, i, j, params, is.data.noisy){
  #function to calculate the emission probability for the velocity jump HMM model
  #assume a two state model, so dependence on current and previous states only
  #calls cpp code for faster calculation of emission probabilities

  #angle_change is a single point from the angle change time series
  #previous_angle_change: obersved angle change from preceeding time interval
  #i: previous hidden state (transition or no transition)
  #j: current hidden state 
  #params: contain model params to infer
  #is.data.noisy: logical for whether measurement noise is present, for calculation of emission probabilities
##############################

  if (is.data.noisy){
    case1 <- dnorm( angle_change, mean = 0, sd = params["sigma"],log=TRUE )
    case2 <- q01CPP(angle_change,params['sigma'],1) 
    case3 <- q010CPP(angle_change,previous_angle_change,params['sigma'],1)
    case4 <- q011CPP(angle_change,previous_angle_change,params['sigma'],1) 
    beta <- ifelse(i>0,ifelse(j>0,case4,case3),ifelse(j>0,case2,case1))
  } else {
    case1 <- ifelse(abs(angle_change)>0,-Inf,0) 
    case2 <- p01CPP(angle_change, TRUE) 
    case3 <- p010CPP(angle_change, previous_angle_change, TRUE) 
    case4 <- p011CPP(angle_change, previous_angle_change, TRUE) 
    beta <- ifelse(i>0,ifelse(j>0,case4,case3),ifelse(j>0,case2,case1))
  }

  return(beta)
}
