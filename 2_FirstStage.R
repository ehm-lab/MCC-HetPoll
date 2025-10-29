#####################################################################
#
#                 MCC HetPoll
#             Part 1: First stage, city level
#
#####################################################################

library(dlnm)
library(splines)
library(dplyr)

#-------------------------------------
# Set up empty objects to save results
#-------------------------------------

# Coefficients
coefall <- matrix(NA, nrow(cities), 1, dimnames = list(cities$city))
redall <- vector("list", nrow(cities))

# Vcov matrices
vcovall <- vector("list", nrow(cities))

# Convergence indicator
conv <- rep(NA, nrow(cities))

# To store total PM
mean_pm <- rep(NA, nrow(cities))

# Names
names(vcovall) <- names(conv) <- names(mean_pm) <- cities$city

#-------------------------------------
# Model parameters
#-------------------------------------

maxlagp <- 1
arglagp <- list(fun = "strata") # Equivalent to MA
argvarp <- list(fun = "lin")

cen <- 0

timedf <- 7

#-------------------------------------
# Loop on cities
#-------------------------------------

for(i in seq(length(dlist))) {
  cat(i,"")
  citydat <- dlist[[i]]
  
  # Construct crossbasis for temperature confounding
  cbt <- crossbasis(citydat$tmean, lag = 3, 
    arglag = list(fun = "strata"),
    argvar = list(fun = "bs", 
      knots = quantile(citydat$tmean, c(.1, .75, .9), na.rm = T))
  )
  
  # Construct crossbasis for PM2.5
  cbp <- crossbasis(citydat$pm25, lag = maxlagp, 
    argvar = argvarp, arglag = arglagp) 
  
  # Estimate the model
  model <- glm(death ~ cbp + cbt + dow + 
      ns(date, df = round(timedf * length(date) / 365.25)), 
    data = citydat, family = quasipoisson, na.action = "na.exclude")
  
  # Succesful estimation?
  conv[i] <- model$converged
  
  # Store results for 2nd stage meta-analysis
  redall[[i]] <- crosspred(cbp, model, cen = cen, at = 10) # Overall
  coefall[i,] <- redall[[i]]$allfit
  vcovall[[i]] <- redall[[i]]$allse^2
  mean_pm[i] <- mean(citydat$pm25, na.rm = T)
}

#-------------------------------------
# Save results
#-------------------------------------

# Put 1st stage results into a nice table
cities <- mutate(cities, coef = coefall, var = unlist(vcovall), conv = conv,
  mean_pm = mean_pm)

# Export it
write.csv(cities, file = "Data/cities_dat.csv", quote = F, row.names = F)

# Also export composition data
tot_spec <- do.call(rbind, dlist_spec)
write.csv(tot_spec, file = "Data/PMcomponents.csv", quote = F, row.names = F)