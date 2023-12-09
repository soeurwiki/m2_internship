#   name : splsda.R
#
#   Author (2023)  Safiya ATIA

library("mixOmics")

X <- readRDS("X_plsda.rds")
Y <- readRDS("Y_plsda.rds")

metabo.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.metabo <- perf(metabo.splsda, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
pdf("perf_splsda_metabo.pdf") 
plot(perf.splsda.metabo, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
dev.off() 

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.metabo <- tune.splsda(X, Y, ncomp = 4, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 10, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 cpus = 2) # allow for paralleliation to decrease runtime

pdf("tune_splsda_metabo.pdf") 
plot(tune.splsda.metabo, col = color.jet(4)) # plot output of variable number tuning
dev.off() 

optimal.ncomp <- tune.splsda.metabo$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.metabo$choice.keepX[1:optimal.ncomp]

print(optimal.ncomp)
print(optimal.keepX)

# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)

saveRDS(final.splsda, "final_splsda.rds")