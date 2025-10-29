#####################################################################
#
#                 MCC HetPoll
#             Part 2: Second stage meta analysis
#              experiment with mixtures
#
#####################################################################

library(mixmeta)
library(compositions)
library(zCompositions)
library(fields)
library(maps)
library(mapdata)

#-------------------------------------
#   Load data
#-------------------------------------

# Results from first-stage and city-level data
cities <- read.csv(file = "Data/cities_dat.csv")

# PM composition data
tot_spec <- read.csv(file = "Data/PMcomponents.csv")

#-------------------------------------
#   Preparing constituent objects
#-------------------------------------

# Components and labels
spec_inds <- grep("PM25", colnames(tot_spec))
spec_names <- c("SO4", "NH4", "NO3", "BC", "OC", "SS", "DUST")
spec_labs <- c(expression(SO[4]^{"2-"}), expression(NH[4]^{"+"}), 
  expression(NO[3]^{"-"}), "BC", "OC", "Sea salt", "Dust")
spec_pal <- c(2, 6, 4, 1, 3, 5, 7)

# Mean per city
mean_comp <- tapply(tot_spec, tot_spec$city, function(x){
  # Zero value imputation as in Mart?n-Fern?ndez et al. (2003)
  imp <- multRepl(x[,spec_inds], label = 0, dl = rep(1e-5, 7),
    z.warning = 1, z.delete = F)
  # Transformation to compositional object
  xc <- acomp(imp)
  # Compositional mean (pass through the irl transformation)
  mean(xc)
})
mean_comp <- do.call(rbind, mean_comp)
p <- ncol(mean_comp)
colnames(mean_comp) <- spec_names

# Reordering
mean_comp <- mean_comp[cities$city,]

# Transforming as ALR
alr_comp <- alr(mean_comp)

#-------------------------------------
#   PC Scores of MCC indicators
#-------------------------------------

# Get indicators
indic_names <- c("oldpopprop", "GDP", "avgtmean", 
  "totalrange", "E_GR_AV00", "E_GR_AV14", "B00", "B15")
indic_form <- sprintf("~ %s", paste(indic_names, collapse = " + "))

# PCA
pca_indic <- prcomp(as.formula(indic_form), data = cities,
  na.action = na.exclude, scale. = T)

# The first two components account for 58% of the variance
indicators <- pca_indic$x[,1:2]

#-------------------------------------
#  Main meta-regression model
#-------------------------------------

# Model with ALR as meta-predictors
metamod <- mixmeta(coef ~ alr_comp + indicators, var, 
  random = ~ 1|country/city, data = cities, method = "reml", subset = conv)

# Q and I2 statistics  
summary(metamod) 

#-------------------------------------
# BLUPS
#-------------------------------------

# Compute RRs
BLUPall <- rep(NA, nrow(cities))
BLUPall[cities$conv] <- exp(blup(metamod))

# Extremes
cities[which.max(BLUPall),]
cities[which.min(BLUPall),]

# Number of RR > 1
sum(BLUPall > 1, na.rm = T)

#-------------------------------------
#  Get meta-coefficients
#-------------------------------------

# Coefficient retrieving
coef_est <- rep(NA, p)
coef_est[-p] <- coef(metamod)[2:p]
coef_est[p] <- -sum(coef_est, na.rm = T) # beta7 = - sum(beta1, ..., beta6)

# Standard errors
se_est <- rep(NA, p)
se_est[-p] <- sqrt(diag(vcov(metamod))[2:p])
se_est[p] <- sqrt(sum(vcov(metamod)[2:p,2:p]))

# Confidence intervals
lo_est <- coef_est - 1.96 * se_est
up_est <- coef_est + 1.96 * se_est
sig_est <- lo_est > 0 | up_est < 0

# We mulitply everything by ln(2) to interpret it as the impact of doubling 
#   the relative proportion of the components
coef_est <- log(2) * coef_est
lo_est <- log(2) * lo_est
up_est <- log(2) * up_est

#-------------------------------------
#  Prediction of the RR for each components separately
#-------------------------------------

# We create a sequence for the component
cseq <- seq(.01, .99, by = .01)

# Overall mean of composition
ov_mean <- mean(acomp(mean_comp))

# Prepare prediction data.frame
newdat <- list(indicators = matrix(0, length(cseq), 2, 
  dimnames = list(NULL, sprintf("PC%i", 1:2))))

# Prepare objects to store predictions, confidence limits and 
#   whether the proportion value is observed for the component
preds <- plo <- pup <- is_obs <- matrix(NA, length(cseq), p)
for (j in seq_len(p)){
  # Create a compositional grid. The component of interests varies 
  #   from 0 to 1 (excluded) and the other are taken as the overall 
  #   mean adjusted for closure while keeping the subcomposition constant
  xj <- matrix(NA, length(cseq), p, dimnames = list(NULL, spec_names))
  xj[,-j] <- sapply(ov_mean[-j] / sum(ov_mean[-j]), "*", 1 - cseq)
  xj[,j] <- cseq
  
  # Transform by ALR and add to the newdat data.frame
  newdat$alr_comp <- alr(xj)
  
  # Make predictions
  pred_obj <- predict(metamod, newdata = newdat, ci = T)
  
  # Store predictions and prediction intervals only in the range of 
  #   observed values
  obs_rng <- range(mean_comp[,j])
  is_obs[,j] <- cseq >= obs_rng[1] & cseq <= obs_rng[2]
  
  # Exponentiate predictions to put on RR scale
  preds[,j] <- exp(pred_obj[,1])
  plo[,j] <- exp(pred_obj[,2])
  pup[,j] <- exp(pred_obj[,3])
}

# Copy preds and keep only observed proportions
preds_obs <- preds
preds_obs[!is_obs] <- NA

#-------------------------------------
#  Comparison with nested models
#-------------------------------------

# Apply model with only indicator PC
metaindic <- mixmeta(coef ~ indicators, var, random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)

# Apply null model: no meta-predictor
metanull <- mixmeta(coef, var, random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)

# Put together
allmodels <- list(full = metamod, indic = metaindic, null = metanull)

#--- create comparison table ---
compar_tab <- data.frame(
  qstat = sapply(allmodels, function(x) summary(x)$qstat$Q),
  i2stat = sapply(allmodels, function(x) summary(x)$i2stat)
)

fwald <- function(full, null) {
  ind <- !names(coef(full)) %in% names(coef(null))
  coef <- coef(full)[ind]
  vcov <- vcov(full)[ind,ind]
  waldstat <- coef %*% solve(vcov) %*% coef
  df <- length(coef)
  pval <- 1 - pchisq(waldstat, df)
  return(list(waldstat = waldstat, pvalue = pval))
}

full_wald <- fwald(metamod, metaindic)
indic_wald <- fwald(metaindic, metanull)

# Add to the table
compar_tab$Wald_stat <- c(full_wald$waldstat, indic_wald$waldstat, NA)
compar_tab$Wald_pvalue <- c(full_wald$pvalue, indic_wald$pvalue, NA)

#---- Export Table 2
write.table(compar_tab, file = "Results/Table2.csv", quote = F,
  row.names = F, sep = ",")

