
library(qgraph)
library(mgm)
library(RColorBrewer)
library(scales)

# All code are adapted from:
# Haslbeck, J. M., Bringmann, L. F., & Waldorp, L. J. (2021). 
# A tutorial on estimating time-varying vector autoregressive models.
# Multivariate behavioral research, 56(1), 120-149.
# https://github.com/jmbh/tvvar_paper

theta <- 0.6 # Maximum parameter value
N <- 1000 # Number of time points
p <- 4 # Number of nodes
x <- seq(from=N*0.2, to=N, length=N*0.8) # Input to Sigmoid function


# Time-varying edge types: constant effect (i.e., stationary) until hypothetical PA intervention
# at time = 200 / 1000, then time-varying based on Sigmoid function
edgetypes <- list()
# Edge 1: Constant effect of -0.3 throughout all time points
edgetypes[[1]] <- rep(-theta/2, N) 
# Edge 2: Constant -0.6 effect until t=200, then linear increase to 0.6
edgetypes[[2]] <- c(rep(-theta, N*0.2),
                    seq(from=-theta, to=theta, length=N*0.8))
# Edge 3: Constant 0.6 effect until t=200, then Sigmoidal decrease
edgetypes[[3]] <- c(rep(max(theta / (1 + exp (0.01 * (x - 500)))), N*0.2),
                    theta / (1 + exp (0.01 * (x - 500))))

# Plot
names(edgetypes) <-c("Constant", "Linear increase", "Sigmoidal decrease")

par(mfrow=c(2, 3))
for(i in 1:3){
  plot(edgetypes[[i]], type="l", lwd="2", ylim=c(-0.6, .6),
       main=names(edgetypes)[i], 
       xlab="Time point", ylab="Parameter Value")
  abline(h = 0, col = "grey", lty = "dashed")
}



# Variable names
labels <- c("PA", "pos_aff", "stress", "neg_aff")

# Matrix with 1's as indicators for non-zero parameter values
ind_mat <- matrix(c(0, 1, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 0, 1,
                    1, 0, 1, 1), nrow=4, byrow=T)

pdf(file="supp_2_plots/qgraph_parameters.pdf")
qgraph(ind_mat, labels=labels, layout="circle", theme="colorblind")
dev.off()

# Table containing location in 4x4 matrix for each effect
effect_ind <- which(ind_mat == 1, arr.ind=T) |> data.frame()

effect_ind$edgetype <- c(1, 2, 3, 3, 3)

effect_ind$tv_param <- c("Constant effect",
                         "Linear increase",
                         "Sigmoidal decrease",
                         "Sigmoidal decrease",
                         "Sigmoidal decrease")

effect_ind$effect <- c("neg_aff -> PA",
                       "PA -> pos_aff",
                       "neg_aff -> stress",
                       "stress -> neg_aff",
                       "neg_aff -> neg_aff")


knitr::kable(effect_ind)

# Create 3 dimensional array (4x4x250) containing directed (autoregressive and cross-lagged) effects of each variable 
# and parameter values of those effects varying across N time points.
B <- array(0, dim=c(p, p, N))

for(i in 1:nrow(effect_ind)){
  B[effect_ind[i, 1], effect_ind[i, 2], ] <- edgetypes[[effect_ind[i, 3]]]
}

# Plot all 4x4 parameters values over time
# pdf(file="supp_2_plots/tv_param_plots.pdf")
# par(mfrow=c(4, 4))
# for(i in 1:p){
#   for(j in 1:p){
#     plot(B[i, j, ],
#          main=paste0(labels[i], " -> ", labels[j]),
#          col="darkred",
#          type="l", lwd="1", ylim=c(-0.7, 0.7),
#          xlab="Time point", ylab="Parameter value"
#     )
#     # Dashed line at y=0
#     abline(h = 0, col = "grey", lty = "dashed")
#   }
# }
# dev.off()

# Plot effects of interest (neg_aff -> PA, PA -> pos_aff, stress -> neg_aff)

titles <- c(expression("(A) Neg_aff"["t-1"]  %->%  "PA"["t"]),
            expression("(B) PA"["t-1"]  %->%  "Pos_aff"["t"]),
            NA,
            expression("(C) Stress"["t-1"]  %->%  "Neg_aff"["t"]))

pdf(file="Fig_3.pdf", height=3, width=9)

lmat <- matrix(c(1, 2, 3), ncol=3, byrow = T)
layout(lmat, heights = c(1, 1, 1), 
       widths = c(1, 1, 1), respect = T)

for(i in c(1, 2, 4)){
  plot(B[effect_ind[i, 1], effect_ind[i, 2], ],
       main=titles[i],
       col="darkred",
       type="l", lwd="1", ylim=c(-0.7, 0.7),
       xlab="Time point", ylab="Parameter value"
    )
    # Dashed line at y=0
    abline(h = 0, col = "grey", lty = "dashed")
}
dev.off()






# Generate data from time-varying VAR
d <- matrix(NA, N, p)

set.seed(1)
d[1, ] <- rnorm(p, mean=0, sd=0.2) # intialized first row

for(t in 2:N){
  for(node in 1:p){
    d[t, node] <- sum(d[t-1, ] * B[node, , t]) + rnorm(1, mean=0, sd=0.2)
  }
}


# Select appropriate bandwidth value using cross-validation to minimize prediction error
bw_seq <- seq(0.01, 1, length = 10)
set.seed(1)
bw_sel <- bwSelect(data = d,
                   type = rep("g", 4),
                   level = rep(1, 4),
                   bwSeq = bw_seq,
                   bwFolds = 1,
                   bwFoldsize = 20,
                   modeltype = "mvar",
                   lags = 1,
                   scale = T,
                   timepoints=seq(0, 1, length=N))

bw_seq[which.min(bw_sel$meanError)] # 0.23


# Fit time-varying VAR 
set.seed(1)
tvvar_mod <- tvmvar(data = d,
                    type = rep("g", p),
                    level = rep(1, p),
                    lambdaSel = "CV",
                    timepoints = seq(0, 1, length=N), # equal spacing between time points
                    estpoints = seq(0, 1, length = 20), # 20 estimation points
                    bandwidth = bw_seq[which.min(bw_sel$meanError)],
                    lags = 1,
                    scale = T)

# save(tvvar_mod, file="tvvar_mod.RData")
load("tvvar_mod.RData")



# 95% bootstrapped confidence intervals
rel <- resample(object=tvvar_mod, 
                data=d, 
                nB=20, 
                blocks=10,
                seeds=1:20, 
                quantiles=c(.05, .95))

# save(rel, file="rel.RData")
load("rel.RData")




# Plot time-varying (modulus of) edge weights from the "wadj"
pdf(file="supp_2_plots/tv_wadj_plots.pdf")
par(mfrow=c(4, 4))
for(i in 1:p){
  for(j in 1:p){
    plot(tvvar_mod$wadj[i, j, 1, ],
         main=paste0(labels[i], " -> ", labels[j]),
         type="p", lwd="1", ylim=c(-0.3, 0.8),
         xlab="Time point", ylab="wadj value"
    )
    # Dashed line at y=0
    abline(h = 0, col = "grey", lty = "dashed")
  }
}
dev.off()



# Plot Networks at different estimation points
tvvar_mod$edgecolor[, , , ][tvvar_mod$edgecolor[, , , ] == "darkgreen"] <- c("blue") # colorblind theme
ep_sel <- c(4, 12, 20) # plot parameters at these estimation points




pdf(file = "Fig_4.pdf", height = 8, width = 8)

lmat <- matrix(c(0, 1, 2, 3,
                 4, 4, 4, 4,
                 5, 5, 5, 5), ncol=4, byrow = T)
layout(lmat, 
       heights = c(2.5, 0.5, 3), 
       widths = c(1, 2, 2, 2), respect = T)


for(tp in ep_sel) {
  qgraph(tvvar_mod$wadj[, , 1, tp], 
         layout = "circle",
         edge.color = tvvar_mod$edgecolor[, , 1, tp], 
         labels = labels, 
         maximum = 0.8,
         vsize = 16, 
         esize = 14,
         asize = 14, 
         minimum = 0.1,
         mar = c(2,5,2,5))
}


par(mar=c(2,6.5,0,4))
plot.new()
plot.window(xlim=c(1,20), ylim=c(0, 50))
axis(1, c(1, 4, 8, 12, 16, 20), lwd=2, labels=T)
mtext("Estimation points", at=9.5, side=1, line=-2, cex=0.8)

arrows(4.5, 50, 4, 10, length=0.1, lwd=2)
arrows(11.65, 50, 12, 10, length=0.1, lwd=2)
arrows(19, 50, 20, 10, length=0.1, lwd=2)



# Line plots of time-varying parameters
par_ests <- tvvar_mod$wadj
ind_negative <- which(tvvar_mod$signs == -1, arr.ind = T)
par_ests[ind_negative] <- par_ests[ind_negative] * -1

plot.new()
par(mar = c(0,4,0,4))
plot.window(xlim=c(0, 20), ylim=c(-.5, .75))
axis(2, c(-0.5, -.25, 0, .25, .5,  0.75), lwd=2, las=2)
abline(h = 0, col = "grey", lty=2)
abline(v = 4, col = "grey", lty=1)
title(ylab = "Parameter estimate", cex.lab = 1.2)
text("Hypothetical", x=2, y=.2, cex=1.2, col="red")
text("PA intervention", x=2, y=.15, cex=1.2, col="red")
arrows(3, .1, 3.8, 0.05, length=0.1, lwd=1.5, col="red")

cols <- brewer.pal(5, "Set1")[3:5] # red, blue, green, purple, orange

# Select which parameters to plot
param_plot <- effect_ind[c(1,2,4), 1:2] |> as.matrix()

for(i in 1:nrow(param_plot)) {
  par_row <- param_plot[i, ]
  ## Plot point estimates
  P1_pointest <- par_ests[par_row[1], par_row[2], 1, ]
  lines(1:20, P1_pointest, col = cols[i], lwd = 2, lty=i) 
  ## Plot uncertainty estimates [new shading]
  # Compute CIs
  CIs <- apply(rel$bootParameters[par_row[1], par_row[2], 1, , ], 1, function(x) {
    quantile(x, probs = c(.05, .95))
  } )
  # Plot shading
  polygon(x = c(1:20, 20:1), y = c(CIs[1,], rev(CIs[2,])), col=alpha(colour = cols[i], alpha = .3), border=FALSE)
  # Legend
   legend_labels <- c(expression("Neg_aff"["t-1"]  %->%  "PA"["t"]),
                     expression("PA"["t-1"]  %->%  "Pos_aff"["t"]),
                     expression("Stress"["t-1"]  %->%  "Neg_aff"["t"]))
  
  legend(4.5, .3, 
         legend_labels,
         col = cols, 
         lwd = 2, bty = "n", cex=1.4, lty=1:3)
  
}

dev.off()






