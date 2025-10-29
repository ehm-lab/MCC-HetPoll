#####################################################################
#
#                 MCC HetPoll
#             Part 1: Summary statistics fo compositional data
#
#####################################################################

library(corrplot)

#-------------------------------------
# Figure 1: Map with BLUPS
#-------------------------------------

# Compute average PM2.5
mean_pm <- cities$mean_pm

# Create colorscale based on RR
cutoff <- seq(0.975, 1.025, by = 0.005)
labels <- paste0(paste0(cutoff[-length(cutoff)], "-", cutoff[-1]))
citycat <- cut(BLUPall, cutoff, labels = labels, include.lowest = T)
pal <- tim.colors(length(labels))

# Create point size based on mean PM2.5
ptsiz_rng <- c(.5, 2.5)
mean_scale <- (mean_pm - min(mean_pm)) / (diff(range(mean_pm)))
pt_size <- mean_scale * diff(ptsiz_rng) + ptsiz_rng[1]
pm_scale <- pretty(mean_pm)
size_scale <- (pm_scale - min(mean_pm)) / (diff(range(mean_pm))) * 
  diff(ptsiz_rng) + ptsiz_rng[1]

#---- Draw map
x11(width = 10)
map("worldHires", mar=c(0,0,0,0), col = grey(0.95),
  myborder = 0, fill = T, border = grey(0.5), lwd = 0.3)
# Add points with size = RR and colorscale for mean PM2.5
points(cities$long, cities$lat, pch = 21, cex = pt_size, bg = pal[citycat])
# Scale and legend
map.scale(-5, -50, ratio = F, cex = 0.7, relwidth = 0.1)
rect(par("usr")[1], 45, -129, -70, border = NA, col = "white")
lg <- legend(-175, 40, labels, pt.cex = 1.2, bg = "white",
  pch = 21, pt.bg = pal, box.col = "white", cex = 0.7, inset = 0.02,
  title = expression(paste("Predicted RR (10 ", mu, "g/", m^3, ")")), 
  title.adj = 0
)
legend(lg$rect$left, lg$rect$top - lg$rect$h, pm_scale, inset = 0.02, 
  pch = 21, pt.cex = size_scale, pt.bg = "grey",box.col = "white", cex = 0.7,   
  title = expression(paste("Mean ", PM[2.5], " concentration (", 
    mu, "g/", m^3, ")")),
  bg = "white", title.adj = 0, y.intersp = 1.5, xjust = 0, 
  ncol = length(pm_scale)
)

dev.print(pdf, file = "Results/Figure1.pdf")

#------------------------------------------
# Figure 2: Mean composition per year and country
#------------------------------------------

# Proportion of zeros of each component
prop_zero <- apply(tot_spec[,spec_inds], 2, function(x) mean(x == 0))

# Zeros imputation
imp_spec <- multRepl(tot_spec[,spec_inds], label = 0, dl = rep(1e-5, 7),
  z.warning = 1, z.delete = F)

# Split by country
count_split <- split(as.data.frame(imp_spec), 
  list(country = tot_spec$country, year = tot_spec$year))
countries <- unique(cities[, c("country", "countryname")])

# Average by country and year
agg_spec <- t(sapply(count_split, function(x) mean(acomp(x))))
agg_spec_rows <- Reduce(rbind,strsplit(rownames(agg_spec), "\\."))

#----- Plot the average
x11(width = 15, height = 10)
par(mfrow = n2mfrow(nrow(countries) + 1, asp = 1.5), 
  mar = c(3, 2, 3, 1))
for (i in seq_len(nrow(countries))){
  # Select country
  country_inds <- agg_spec_rows[,1] == countries$country[i]
  country_dat <- agg_spec[country_inds,]
  country_dat <- merge(cbind(as.numeric(agg_spec_rows[country_inds, 2]), country_dat), 
    2003:2017, by = 1, all.y = T)
  
  # Barplot
  bp <- barplot(t(data.matrix(country_dat[,-1])), col = spec_pal, 
    border = NA, names.arg = country_dat[,1], 
    main = countries$countryname[i])
}
plot.new()
legend("topleft", spec_labs, fill = spec_pal, bty = "n", ncol = 3,
  cex = 1.8, xpd = NA, border = NA)

dev.print(pdf, file = "Results/Figure2.pdf")

#-------------------------------------
#  Figure 3: Meta-regression resuls
#-------------------------------------

# Panel matrix
pm <- rbind(c(9, rep(10, 3)), 
  cbind(1, matrix(c(2:7, 0, 8, 0), nrow = 3, byrow = T)))

# Initialize plot
x11(height = 10, width = 15)
par(mar = c(5, 6, 3, 1), cex.main = 1.5, cex.lab = 1.2)
layout(pm, width = c(.4, .2, .2, .2), height = c(.1, .3, .3, .3))

#----- Panel A: meta-coefficients

plot(exp(coef_est), -seq_len(7), 
  xlab = "Relative Excess Risk", ylab = "", cex.lab = 1.3,
  axes = F, xlim = range(c(exp(lo_est), exp(up_est))))
abline(v = 1, lty = 2)
segments(exp(lo_est), -seq_len(7), exp(up_est), -seq_len(7), lwd = 2,
  col = grey(.2))
points(exp(coef_est), -seq_len(7), cex = 3, 
  pch = ifelse(sig_est, 15, 16), col = spec_pal)
axis(1, cex.axis = 1.1)
axis(2, at = -seq_len(7), labels = spec_labs, 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0, cex.axis = 1.5)
box()

#----- Panel B: predicted values
par(mar = c(5, 4, 3, 1))
for (j in seq_len(p)){
  # Initiate an empty plot with labels
  plot(0, 0, col = NA, xlim = 100 * range(mean_comp), 
    ylim = c(min(plo, na.rm = T), max(pup, na.rm = T)),
    ylab = "RR", xlab = "Proportion (%)",
    main = bquote(bold(.(spec_labs[j][[1]]))))
  # Draw prediction interval
  polygon(100 * c(cseq[!is_obs[,j]], rev(cseq[!is_obs[,j]])), 
    c(plo[!is_obs[,j],j], rev(pup[!is_obs[,j],j])), 
    col = adjustcolor(spec_pal[j], .3), border = NA, density = 20)
  polygon(100 * c(cseq[is_obs[,j]], rev(cseq[is_obs[,j]])), 
    c(plo[is_obs[,j],j], rev(pup[is_obs[,j],j])), 
    col = adjustcolor(spec_pal[j], .2), border = NA)
  # Add prediction
  lines(100 * cseq, preds[,j], lwd = 1, col = spec_pal[j])
  lines(100 * cseq[is_obs[,j]], preds[is_obs[,j],j], lwd = 3, col = spec_pal[j])
  # Add lines indicating observed range
  abline(v = 100 * range(cseq[is_obs[,j]]), lty = 2)
  abline(h = 1)
}

#----- Add letters for plots
par(mar = rep(0, 4), cex.main = 1.5, cex.lab = 1.2)
plot.new()
text(par("usr")[1], par("usr")[3], "A", cex = 3, adj = c(0, 0), xpd = T)

plot.new()
text(par("usr")[1], par("usr")[3], "B", cex = 3, adj = c(0, 0), xpd = T)

#----- Save
dev.print(pdf, file = "Results/Figure3.pdf")
  

  