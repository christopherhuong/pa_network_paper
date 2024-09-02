
library(huge)
library(qgraph)
library(bootnet)
library(ggplot2)
library(gridExtra)



# glasso regularized GGM networks from cross-sectional data --------------------------

# Load in data
stretch <- read.csv('StretchingGroupData.csv')
cycle <- read.csv('CyclingGroupData.csv')

# Nonparanormal transformation to 
stretch <- huge.npn(stretch, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
stretch <- as.data.frame(apply(stretch, use='pairwise.complete.obs', 2, as.numeric))

cycle <- huge.npn(cycle, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
cycle <- as.data.frame(apply(cycle, use='pairwise.complete.obs', 2, as.numeric))


names = c("Perserv", "Neg", "Critic", "Brood", "Replay")
nodenames =c("Perseveration; “How much have you been thinking about your performance?”",
             "Negativity; “How negative have your thoughts been?",
             "Self-criticism; “How much have you criticized yourself about your performance?",
             "Brooding; “How much have you thought about how upset you felt?”",
             "Replaying; “To what extent did you replay parts of what happened in your mind?")

# Estimate glasso networks 
# Tuning parameter (lambda) of the glasso regularization is set by minimizing the
# EBIC, with its hyperparameter (gamma) set to 0.5
# gamma closer to 1 = higher specificty
# gamma closer to 0 = higher sensitivity
stretch_glasso <- estimateNetwork(stretch, default = "EBICglasso", tuning = 0.5)
cycle_glasso <- estimateNetwork(cycle, default = "EBICglasso", tuning = 0.5)

L <- averageLayout(stretch_glasso, cycle_glasso)

# Save layouts
stretch_plot <- qgraph(stretch_glasso$graph,
                       vsize=14, 
                       esize=14,
                       title = "(A) Stretching",
                       title.cex=1.5,
                       layout=L, labels=names, theme = "colorblind",
                       DoNotPlot=T)


cycle_plot <- qgraph(cycle_glasso$graph, 
                     vsize=14, 
                     esize=14,
                     title = "Cycling",
                     title.cex=1.5,
                     layout=L, labels=names, theme = "colorblind",
                     DoNotPlot=T)


all_cen <- list("Stretching"=stretch_glasso, "Cycling"=cycle_glasso)
cen_plot <- centralityPlot(all_cen, include = "Strength", labels = names, scale = "z-scores") + xlab("z-scores") + ggtitle("(B)")


# Plot Figure 1 A and B
pdf(file="Fig_1A.pdf", width=8, height=4)
par(mfrow=c(1, 2))
plot(stretch_plot)
plot(cycle_plot)
dev.off()

pdf(file="Fig_1B.pdf", width=4, height=4)
plot(cen_plot)
dev.off()






# Personalized VAR-based networks -----------------------------------------
library(graphicalVAR)
library(qgraph)

n_var <- 4
vars <- c("PA", "depr", "conc", "tired")

# Simulate two random networks
set.seed(321)

# Randomly generate a "true" underlying temporal and contemporaneous network structure for each individual
A <- randomGVARmodel(Nvar = n_var, probKappaEdge = 0.2, probKappaPositive = 0.5,
                     probBetaEdge = 0.3, probBetaPositive = 0.3)
B <- randomGVARmodel(Nvar = n_var, probKappaEdge = 0.2, probKappaPositive = 0.5,
                     probBetaEdge = 0.3, probBetaPositive = 0.3)

# Simulate multivariate normal data from the beta matrix of lagged effects, and kappa matrix of contemporaneous effects
A_data <- graphicalVARsim(nTime=100,
                          beta = A$beta, # Beta matrix
                          kappa = A$kappa) # Kappa matrix

B_data <- graphicalVARsim(nTime=100, 
                          beta = B$beta, 
                          kappa = B$kappa)

# Fit a graphical VAR model on each persons data
A_GVAR <- graphicalVAR(A_data, 
                       nLambda=50, 
                       gamma=0.5) # Tuning parameter to control level of sparsity.

B_GVAR <- graphicalVAR(B_data, 
                       nLambda=50, 
                       gamma=0.5)


gvar_list <- list(A_PDC = A_GVAR$PDC,
                  A_PCC = A_GVAR$PCC,
                  B_PDC = B_GVAR$PDC,
                  B_PCC = B_GVAR$PCC)

titles <- c("Individual A: Temporal Network", "Individual B: Temporal Network",
            "Individual A: Contemporaneous Network", "Individual B: Contemporaneous Network")
# Plot the model estimates
# Figure 2 ------------------------------------------------------------------
pdf(file="Fig_2.pdf", width=8, height=8)

lmat <- matrix(c(1, 0, 2, 
                 3, 0, 4), ncol=3, byrow = T)
lo <- layout(lmat, 
             heights = c(1, 1), 
             widths = c(1, 0.1, 1), respect = T)

for(i in 1:4){
  qgraph(gvar_list[[i]], layout="circle", theme="colorblind", labels=vars,
         vsize=16, 
         esize=12,
         asize=12,
         title=titles[i], title.cex=1.5,
         mar=c(7, 7, 7, 7))
}
dev.off()











