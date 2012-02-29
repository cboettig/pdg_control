# Bonvin NCO Control under model uncertainty


# target S
# Estimate S_hat

# Use closure to define a profit function
profit <- profit_harvest()

# initialize system
x_h <- numeric(OptTime) # population dynamics with harvest
h <- numeric(OptTime) # havest level
P <- numeric(OptTime) # Profit 
S_hat <- numeric(OptTime) # Escapement (control strategy) 
x_h[1] <- x0  # initial values
S_hat[1] <- 0.2 # initial guess


# next time step
h[t] <- x_h[t]-S_hat[t] # harvest = stock minus escapement
x_h[t+1] <- zg() * f(x_h[t], h[t], pars) # with harvest

# calculate cost -- huh?  
profit(x_h[t+1], h[t])
