#--------------------#
# CLEAN ENVIRORNMENT #
#--------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
# Option 1)
# You can type the path to the directory 
# where you have saved this R script.
# E.g.:
#
# wd <- "~/Workshop/BayesianInference/R"
# setwd( wd )
#
# Option 2)
# You can find your working directory with the next
# commands. If you are not going to use this
# option and prefer option 1), do not run the
# following commands or just comment them:
filename <- "Practical_1.R"
filepath <- file.choose( )  # Browse and select
#                           # "Practical_1.R"
#                           # in the window
wd       <- substr( filepath, 1, nchar( filepath ) - nchar( filename ) )
# If you are a Windows user, you might want to run
# the next command in case there are issues with 
# the backslash
wd <- gsub( pattern = "[\\]", replace = "/", x = wd )
setwd( wd )
#
# Option 3)
# Alternatively, you can install a package that
# automatically finds the path to this R script.
# You can comment line 38 if you have already
# installed the `rstudioapi` package
install.packages( "rstudioapi" )
library( rstudioapi ) 
# The command in line 43 will get the path to this
# R script and will assign it to the variable
# `path_to_file`
path_to_file <- getActiveDocumentContext()$path 
# The next commands will use the variable
# `path_to_file` to locate the path to the
# working directory
wd      <- paste( dirname( path_to_file ), "/", sep = "" )
wd.name <- dirname( path_to_file )
# Last, we use the command in line 52 to set
# the working directory to the path saved in
# object `wd`
setwd( wd )

########################
# PART 1: Introduction #
########################
# The data are the 12S rRNA alignment of human and orangutan, with 948 base
# pairs and 90 differences (84 transitions and 6 transversions).
# Example data from Table 1.3, p.7 of Yang (2014), Molecular Evolution:
# A Statistical Approach. Oxford University Press.

# 1. Define data variables
# Length of alignment in bp
n  <- 948
# Total number of transitions (23+30+10+21)
ns <- 84
# Total number of transversions (1+0+2+0+2+1+0+0)
nv <- 6

# 2. Define log-likelihood function, f(D|d,k)
#    This function uses Kimura's (1980) substitution model -- see p.8 in Yang (2014)
#
# Arguments:
#
#   d  Numeric, value for the distance.
#   k  Numeric, value for parameter kappa.
#   n  Numeric, length of the alignment. Default: 948.
#   ns Numeric, total number of transitions. Default: 84.
#   nv Numeric, total number of transversions. Default: 6.
k80.lnL <- function( d, k, n = 948, ns = 84, nv = 6 ) {

  # Define probabilities
  p0 <- .25 + .25 * exp( -4*d/( k+2 ) ) + .5 * exp( -2*d*( k+1 )/( k+2 ) )
  p1 <- .25 + .25 * exp( -4*d/( k+2 ) ) - .5 * exp( -2*d*( k+1 )/( k+2 ) )
  p2 <- .25 - .25 * exp( -4*d/( k+2 ) )

  # Return log-likelihood
  return ( ( n - ns - nv ) * log( p0/4 ) +
            ns * log( p1/4 ) + nv * log( p2/4 ) )

}

# 3. Define variables to plot
# Dimension for the plot
dim <- 100
# Vector of d values
d.v <- seq( from = 0, to = 0.3, len = dim )
# Vector of k values
k.v <- seq( from = 0, to = 100, len = dim )
dk  <- expand.grid( d = d.v, k = k.v )

# We can compute the log-likelihood surface, f(D|d,k), for 
# each possible combination of d and k that are saved
# in object dk
lnL <- matrix( k80.lnL( d = dk$d, k = dk$k ), ncol = dim )
# For numerical reasons, we scale the likelihood to be 1
# at the maximum, i.e. we subtract max(lnL)
L <- exp( lnL - max( lnL ) )

# Now, we can also compute the prior surface, f(D)f(k), 
# and the unscaled posterior surface, f(d)f(k)f(D|d,k),
# so we can compare them against the likelihood surface
Pri <- matrix( dgamma( x = dk$d, shape = 2, rate = 20 ) *
               dgamma( x = dk$k, shape = 2, rate = .1 ),
               ncol = dim )
Pos <- Pri * L

# 4. Plot prior, likelihood, and unscaled posterior surfaces
#    We want one row and three columns
par( mfrow = c( 1, 3 ) )
# Prior surface. Note that the `contour()` function creates a contour plot
image( x = d.v, y = k.v, z = -Pri, las = 1, col = heat.colors( 50 ),
       main = "Prior", xlab = "distance, d",
       ylab = "kappa, k", cex.main = 2.0,
       cex.lab = 1.5, cex.axis = 1.5 )
contour( x = d.v, y = k.v, z = Pri, nlev=10, drawlab = FALSE, add = TRUE )
# Likelihood surface + contour plot
image( x = d.v, y = k.v, z = -L, las = 1, col = heat.colors( 50 ),
       main = "Likelihood", xlab = "distance, d",
       ylab = "kappa, k", cex.main = 2.0,
       cex.lab = 1.5, cex.axis = 1.5 )
contour( x = d.v, y = k.v, z = L, nlev = 10,
         drawlab = FALSE, add = TRUE)
# Unscaled posterior surface + contour plot
image( x = d.v, y = k.v, z = -Pos, las = 1, col = heat.colors( 50 ),
       main = "Posterior", xlab = "distance, d",
       ylab = "kappa, k", cex.main = 2.0,
       cex.lab = 1.5, cex.axis = 1.5 )
contour( x = d.v, y = k.v, z = Pos, nlev = 10,
         drawlab = FALSE, add = TRUE )


###########################################
# PART 2: Markov Chain Monte Carlo (MCMC) #
###########################################
# Now, we want to obtain the posterior distribution by MCMC sampling.
# In most practical problems, constant z cannot be calculated (either
# analytically or numerically), and so the MCMC algorithm becomes necessary.

# 1. Define function that returns the logarithm of the unscaled posterior:
#                             f(d) * f(k) * f(D|d,k)
#    By, default we set the priors as:
#                  f(d) = Gamma(d | 2, 20) and f(k) = Gamma(k | 2, .1)
#
# Arguments:
#
#   d     Numeric, value for the distance.
#   k     Numeric, value for parameter kappa.
#   n     Numeric, length of the alignment. Default: 948.
#   ns    Numeric, total number of transitions. Default: 84.
#   nv    Numeric, total number of transversions. Default: 6.
#   a.d.  Numeric, alpha value of the Gamma distribution that works as a prior
#         for the distance (d). Default: 2.
#   b.d.  Numeric, beta value pf the Gamma distribution that works as a prior
#         for parameter distance (d). Default: 20.
#   a.k.  Numeric, alpha value for the Gamma distribution that works as a prior
#         for parameter kappa (k). Default: 2.
#   b.k.  Numeric, beta value for the Gamma distribution that works as a prior
#         for parameter kappa (k). Default: 0.1.
ulnPf <- function( d, k, n = 948, ns = 84, nv = 6,
                   a.d = 2, b.d = 20, a.k = 2, b.k = .1 ){

  # The normalizing constant in the prior densities can be ignored
  lnpriord <- ( a.d - 1 )*log( d ) - b.d * d
  lnpriork <- ( a.k - 1 )*log( k ) - b.k * k

  # Define log-Likelihood (K80 model)
  expd1 <- exp( -4*d/( k+2 ) )
  expd2 <- exp( -2*d*( k+1 )/( k+2 ) )
  p0 <- .25 + .25 * expd1 + .5 * expd2
  p1 <- .25 + .25 * expd1 - .5 * expd2
  p2 <- .25 - .25 * expd1
  lnL <- ( ( n - ns - nv ) * log( p0/4 ) + ns * log( p1/4 ) + nv * log( p2/4 ) )

  # Return unnormalised posterior
  return ( lnpriord + lnpriork + lnL )
}

# 2. Draft function with MCMC algorithm:
#      1. Set initial states for d and k.
#      2. Propose a new state d* (from an appropriate proposal density).
#      3. Accept or reject the proposal with probability:
#            min(1, p(d*)p(x|d*) / p(d)p(x|d)).
#            If the proposal is accepted set d = d*, otherwise d = d.
#      4. Save d.
#      5. Repeat 2-4 for k.
#      6. Go to step 2.
#
# Arguments:
#
#   init.d  Numeric, initial state value for parameter d.
#   init.k  Numeric, initial state value for paramter k.
#   N       Numeric, number of iterations that the MCMC will run.
#   w.d     Numeric, width of the sliding-window proposal for d.
#   w.k     Numeric, width of the sliding-window proposal for k.
mcmcf <- function( init.d, init.k, N, w.d, w.k ) {

  # We keep the visited states (d, k) in sample.d and sample.k
  # for easy plotting. In practical MCMC applications, these
  # are usually written into a file. These two objects are numeric
  # vectors of length N+1.
  sample.d <- sample.k <- numeric( N+1 )

  # STEP 1: Set initial parameter values to be used during the first
  #         iteration of the MCMC
  # 1.1. Get initial values for parameters k and d. Save these values
  #      in vectors sample.d and sample.k
  d <- init.d;  sample.d[1] <- init.d
  k <- init.k;  sample.k[1] <- init.k
  # 1.2. Get unnormalised posterior
  ulnP  <- ulnPf( d = d, k = k )
  # 1.3. Initialise numeric vectors that will be used to keep track of
  #      the number of times proposed values for each parameter,
  #      d and k, have been accepted
  acc.d <- 0; acc.k <- 0
  # 1.4. Start MCMC, which will run for N iterations
  for ( i in 1:N ){

    # STEP 2: Propose a new state d*
    #         We use a uniform sliding window of width w with reflection
    #         to propose new values d* and k*
    # 2.1. Propose d* and accept/reject the proposal
    dprop <- d + runif( n = 1, min = -w.d/2, max = w.d/2 )
    # 2.2. Reflect if dprop is negative
    if ( dprop < 0 ) dprop <- -dprop
    # 2.3. Compute unnormalised posterior
    ulnPprop <- ulnPf( d = dprop, k = k )
    lnalpha  <- ulnPprop - ulnP

    # STEP 3: Accept or reject the proposal:
    #            if ru < alpha accept proposed d*
    #            else reject and stay where we are
    if ( lnalpha > 0 || runif( n = 1 ) < exp( lnalpha ) ){
      d      <- dprop
      ulnP   <- ulnPprop
      acc.d  <- acc.d + 1
    }

    # STEP 4: Repeat steps 2-3 to propose a new state k*
    # 4.1. Propose k* and accept/reject the proposal
    kprop <- k + runif( n = 1, min = -w.k/2, max = w.k/2 )
    # 4.2. Reflect if kprop is negative
    if ( kprop < 0 ) kprop <- -kprop
    # 4.3. Compute unnormalised posterior
    ulnPprop <- ulnPf( d = d, k = kprop )
    lnalpha  <- ulnPprop - ulnP
    # 4.4. Accept/reject proposal:
    #          if ru < alpha accept proposed k*
    #          else reject and stay where we are
    if ( lnalpha > 0 || runif( n = 1 ) < exp( lnalpha ) ){
      k     <- kprop
      ulnP  <- ulnPprop
      acc.k <- acc.k + 1
    }

    # STEP 5: Save chain state for each parameter so we can later
    #         plot the corresponding histograms
    sample.d[i+1] <- d
    sample.k[i+1] <- k
  }

  # Print out the proportion of times
  # the proposals were accepted
  cat( "Acceptance proportions:\n", "d: ", acc.d/N, " | k: ", acc.k/N, "\n" )

  # Return vector of d and k visited during MCMC
  return( list( d = sample.d, k = sample.k ) )

}

# 2. Now, we are going to set a seed number so we can 
# later reproduce the results that we get with the MCMCs
# that are to be run from this part of the tutorial 
# until the end. You can ommit running the next command
# if you want to get different results each time you 
# run this tutorial
set.seed( 12345 )

# 3. Test run-time
system.time( mcmcf( init.d = 0.2, init.k = 20, N = 1e4,
                    w.d = .12, w.k = 180 ) )
# With PC i7-8750H CPU, 16.0GB RAM, 2.20GHz,
# this MCMC takes ~0.2s to finish

# 4. Run again and save MCMC output
dk.mcmc <- mcmcf( init.d = 0.2, init.k = 20, N = 1e4,
                  w.d = .12, w.k = 180 )

############################
# PART 3: MCMC diagnostics #
############################
# [[ TRACE PLOTS ]] #
# 1. Plot traces of the values sampled for each parameter
par( mfrow = c( 1,3 ) )
# Plot trace for parameter d
plot( x = dk.mcmc$d, ty = 'l', xlab = "Iteration",
      ylab = "d", main = "Trace of d" )
# Plot trace for parameter k
plot( x = dk.mcmc$k, ty = 'l', xlab = "Iteration",
      ylab = "k", main = "Trace of k" )
# Plot the joint sample of d and k (points sampled from posterior surface)
plot( x = dk.mcmc$d, y = dk.mcmc$k, pch = '.', xlab = "d",
      ylab = "k", main = "Joint of d and k" )

# [[ AUTOCORRELATION FUNCTION (ACF) PLOTS ]] #
# Values sampled in an MCMC chain are autocorrelated because new states
# are either the previous state or a modification of it.
# The efficiency of an MCMC chain is closely related to the autocorrelation.
# Intuitively, if the autocorrelation is high, the chain will be inefficient,
# i.e., we will need to run the chain for a long time to obtain a good
# approximation to the posterior distribution.
# The efficiency of a chain is defined as:
#                  eff = 1 / (1 + 2(r1 + r2 + r3 + ...))
# where r_{i} is the correlation for lag i.

# 1. Run a very long chain to calculate efficiency
#    NOTE: With PC i7-8750H CPU, 16.0GB RAM, 2.20GHz,
#          the MCMC takes ~13s to finish with the settings specified in
#          the function below
dk.mcmc2 <- mcmcf( init.d = 0.2, init.k = 20, N = 1e6,
                   w.d = .12, w.k = 180 )

# 2. Plot ACF (the autocorrelation function) for each parameter
par( mfrow = c( 1,2 ) )
acf( x = dk.mcmc2$d )
acf( x = dk.mcmc2$k )

# [[ CHAIN EFFICIENCY ]] #
# Apart from calculating and plotting the ACF, we can also 
# calculate chain efficiency. For this purpose, we use the 
# following function:

# 1. Define efficiency function
#
# Arguments:
#  acf  Numeric, autocorrelation value
eff <- function( acf ) 1 / ( 1 + 2 * sum( acf$acf[-1] ) )

# 2. Compute efficiency for each parameter when running the MCMC
#    with the parameter values set above:
#       mcmcf( init.d = 0.2, init.k = 20, N = 1e6, w.d = .12, w.k = 180 )
eff( acf = acf( dk.mcmc2$d ) )
eff( acf = acf( dk.mcmc2$k ) )
# The efficienciy values are roughly 22% and 20% for d and k respectively

# [[ TIME TO PRACTICE ]] #
# Now, to illustrate inefficient chains, we will run our MCMC again by using
# a proposal density with a too large step size for d, and another with a
# too small step size for k if compared to the much more optimal values
# used before:
dk.mcmc3 <- mcmcf( init.d = 0.2, init.k = 20, N = 1e4,
                   w.d = 3, w.k = 5 )

# 1. Plot traces for each parameter. Note that, because proposal width for d
#    is too large, chain gets stuck at same values of d. On the other hand,
#    proposal width for k is too small, so chain moves slowly
par( mfrow = c( 1,2 ) )
plot( x = dk.mcmc3$d, ty = 'l', main = "Trace of d", cex.main = 2.0,
      cex.lab = 1.5, cex.axis = 1.5, ylab = "d" )
plot( x = dk.mcmc3$k, ty = 'l', main = "Trace of k", cex.main = 2.0,
      cex.lab = 1.5, cex.axis = 1.5, ylab = "k" )

# Now, we run the chain longer but keep the same starting values for the
# rest of the parameters. Then, we compute the chain efficiency for each
# parameter:
dk.mcmc4 <- mcmcf( init.d = 0.2, init.k = 20, N = 1e6,
                   w.d = 3, w.k = 5 )
eff( acf = acf( dk.mcmc4$d, lag.max = 2e3 ) )
eff( acf = acf( dk.mcmc4$k, lag.max = 2e3 ) )
# Efficiency values are roughly 1.5% for d and 0.35% for k

# 1. Plot the traces for efficient (part 2) and inefficient chains
par( mfrow = c( 2,2 ) )
plot( dk.mcmc$d, ty = 'l', las = 1, ylim = c( .05,.2 ),
      main = "Trace of d, efficient chain", xlab = '',
      ylab = "Distance, d", cex.main = 2.0, cex.lab = 1.5 )
plot( dk.mcmc3$d, ty = 'l', las = 1, ylim = c( .05,.2 ),
      main = "Trace of d, inefficient chain", xlab='',
      ylab = '', cex.main = 2.0, cex.lab = 1.5 )
plot( dk.mcmc$k, ty = 'l', las = 1, ylim = c( 0,100 ),
      main = "Trace of k, efficient chain",
      xlab = '', ylab = "ts/tv ratio, k",
      cex.main = 2.0, cex.lab = 1.5 )
plot( dk.mcmc3$k, ty = 'l', las = 1, ylim = c( 0,100 ),
      main = "Trace of k, inefficient chain",
      xlab = '', ylab = '', cex.main = 2.0, cex.lab = 1.5 )

# We will now illustrate the concept of burn-in.
# We will run a chain with a high starting value,
# and another with a low starting value for parameters d and k.
dk.mcmc.l <- mcmcf( init.d = 0.01, init.k = 20, N = 1e4,
                    w.d = .12, w.k = 180 )
dk.mcmc.h <- mcmcf( init.d = 0.4, init.k = 20, N = 1e4,
                    w.d = .12, w.k = 180 )

# 1. Compute man and sd of d.
# We use the low chain to calculate the mean and sd of d, but we
# could also have used the high chain
mean.d <- mean( dk.mcmc.l$d )
sd.d   <- sd( dk.mcmc.l$d )

# 2. Plot the chains
plot( dk.mcmc.l$d, xlim = c( 1,200 ), ylim = c( 0,0.4 ), ty = "l" )
lines( dk.mcmc.h$d, col = "red" )
# Plot a horizontal dashed line to indicate (approximately)
# the 95% CI.
abline( h = mean.d + 2 * c( -sd.d, sd.d ), lty = 2 )
# NOTE: You can use this plot to see how the chains move from either
# the high or low starting values towards the stationary phase (the area
# within the dashed lines). The area before it reaches stationarity is the
# burn-in.

# Last, let's compare chain efficiency.
# Run an efficient chain (i.e., good proposal step sizes)
dk.mcmc.b <- mcmcf( init.d = 0.05, init.k = 5, N = 1e4,
                    w.d = .12, w.k = 180 )
# Run an inefficient chain (i.e., bad proposal step sizes)
dk.mcmc3.b <- mcmcf( init.d  = 0.05, init.k = 5, N = 1e4,
                     w.d = 3, w.k = 5 )

# 1. Plot and compare histograms
# Set breaking points for the plot
bks <- seq( from = 0, to = 150, by = 1 )
# Start plotting
par( mfrow = c( 1,2 ) )
hist( x = dk.mcmc.b$k, prob = TRUE, breaks = bks, border = NA,
      col = rgb( 0, 0, 1, .5 ), las = 1, xlab = "kappa",
      xlim = c( 0,100 ), ylim = c( 0,.055 ) )
hist( x = dk.mcmc$k, prob=TRUE, breaks=bks, border=NA,
      col=rgb(.5, .5, .5, .5), add=TRUE)
hist( x = dk.mcmc3.b$k, prob=TRUE, breaks=bks, border=NA,
      col=rgb(0, 0, 1, .5), las=1, xlab="kappa",
      xlim=c(0,100), ylim=c(0,.055))
hist( x = dk.mcmc3$k, prob=TRUE, breaks=bks, border=NA,
      col=rgb(.5, .5, .5, .5), add=TRUE)

# 2. Calculate the posterior means and s.d for each chain
# Compute means for efficient chains (they are quite similar)
mean( dk.mcmc$d ); mean( dk.mcmc.b$d )
mean( dk.mcmc$k ); mean( dk.mcmc.b$k )
# Compute means for inefficient chains (not so similar)
mean( dk.mcmc3$d ); mean( dk.mcmc3.b$d )
mean( dk.mcmc3$k ); mean( dk.mcmc3.b$k )
# Standard error of the means for efficient chains
sqrt( 1/1e4 * var( dk.mcmc$d ) / 0.23 ) # roughly 2.5e-4
sqrt( 1/1e4 * var( dk.mcmc$k ) / 0.20 ) # roughly 0.2
# Standard error of the means for inefficient chain
sqrt( 1/1e4 * var( dk.mcmc3$d ) / 0.015 ) # roughly 9.7e-4
sqrt( 1/1e4 * var( dk.mcmc3$k ) / 0.003 ) # roughly 1.6

# 3. Plot densities (smoothed histograms) for the efficient and
#    inefficient chains
# Set value to scale the kernel densities for the MCMCs
adj <- 1.5
par( mfrow = c( 1,2 ) )
# Efficient chains
plot( x = density( x = dk.mcmc.b$k, adjust = adj ), col = "blue", las = 1,
      xlim  = c( 0, 100 ), ylim = c( 0, .05 ), xaxs = "i", yaxs = "i" )
lines( x = density( x = dk.mcmc$k, adjust = adj ), col = "black" )
# Inefficient chains
plot( x = density( dk.mcmc3.b$k, adjust = adj ), col = "blue", las = 1,
      xlim = c(0, 100), ylim = c( 0, .05 ), xaxs = "i", yaxs = "i" )
lines( x = density( x = dk.mcmc3$k, adjust = adj ), col = "black" )
