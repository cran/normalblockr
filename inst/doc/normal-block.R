## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----set-up-------------------------------------------------------------------
library(normalblockr)
library(pheatmap)
library(paletteer)

## ----set-seed-----------------------------------------------------------------
set.seed(1)

## ----fix-simulations-parameters-----------------------------------------------
n = 100 # Number of samples (number of rows in the matrix of observations)
p = 40  # Number of entities observed (number of columns in the matrix of observations)
d = 2   # Number of covariates
q = 3   # Number of clusters
kappa = 0 # Mean zero-inflation probability (can also be a vector to define one ZI-probability for each variable). kappa = 0 means that there is no zero-inflation.
omega_structure = "erdos-renyi" # Network structure.
u_v = c(0.3, 0.1) # Parameters to generate an association matrix from a graph, details given in the bibliography.
SNR = 0.75 # Signal to Noise Ratio, defines the relative weight of the covariates and the variance
alpha = rep(1/q, q) # Vector giving probabilities of belonging to each cluster
range_X = c(0, 10)  # Min and max values for the covariates
range_D = c(0.5, 1.5) # Min and max values for the individual entities variances

## ----simulation1--------------------------------------------------------------
my_nb_data <- generate_normal_block_var_data(n, p, d, q, kappa, omega_structure, u_v,
                                            SNR, alpha, range_X, range_D)

## ----simulation1-visualization------------------------------------------------
pheatmap::pheatmap(my_nb_data$Y,
                   color =paletteer::paletteer_c("ggthemes::Orange-Gold", n = 100), cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE)

## ----simulation2--------------------------------------------------------------
kappa_zi <- rnorm(p, mean = 0.7, sd = 0.05)
kappa_zi <- unlist(lapply(kappa_zi, f <- function(x) return(max(0, min(x, 0.9)))))
my_nb_data_zi <- generate_normal_block_var_data(n = n, p = p, d = d, q = 5, kappa = kappa_zi,
                                                omega_structure = omega_structure, u_v = u_v,
                                                SNR = SNR, alpha = alpha,
                                                range_X = range_X, range_D = range_D)

## ----simulation2-visualization------------------------------------------------
min_val <- min(my_nb_data_zi$Y) ; max_val <- max(my_nb_data_zi$Y)
orange_gold_pal <- paletteer::paletteer_c("ggthemes::Orange-Gold", n = 100)
n_breaks <- 100
zero_pos <- round((0 - min_val) / (max_val - min_val) * n_breaks) + 1
custom_pal <- c(orange_gold_pal[1:(zero_pos - 1)], "white", orange_gold_pal[zero_pos:n_breaks])
pheatmap::pheatmap(my_nb_data_zi$Y,
                   color = custom_pal,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = FALSE)

## ----NormalBlockData----------------------------------------------------------
my_data    <- NormalBlockData$new(my_nb_data$Y, my_nb_data$X)
my_data_zi <- NormalBlockData$new(my_nb_data_zi$Y, my_nb_data_zi$X)

## ----NormalBlockData-alt------------------------------------------------------
colnames(my_data$X) <- c("X1", "X2")
my_data_alt    <- NormalBlockData$new(my_data$Y, my_data$X, formula = ~ 0 + X1)

## ----simulation1-NB1----------------------------------------------------------
my_NB <- normal_block(data = my_data,
                      blocks = my_nb_data$parameters$C)

## ----simulation1-NB1-plot-----------------------------------------------------
plot(my_NB)

## ----simulation1-NB1-print----------------------------------------------------
print(my_NB)

## ----simulation1-NB1-plot_network---------------------------------------------
 my_NB$plot_network()

## ----simulation1-NB2----------------------------------------------------------
my_NB <- normal_block(data = my_data,
                      blocks = 3)
print(my_NB)
plot(my_NB)

## ----simulation1-ARI----------------------------------------------------------
aricode::ARI(my_NB$clustering, apply(my_nb_data$parameters$C, 1, which.max))

## ----simulation1-NB3----------------------------------------------------------
my_NB_unknown <- normal_block(data = my_data,
                              blocks = 2:5)

## ----simulation1-NB3-plot-----------------------------------------------------
plot(my_NB_unknown)

## ----simulation1-NB3-model-selection------------------------------------------
myNB_3   <- my_NB_unknown$get_model(3)
myNB_BIC <- my_NB_unknown$get_best_model("BIC")

## ----simulation1-NB-fixed-sparsity--------------------------------------------
my_NB_sparse_low <- normal_block(data = my_data,
                                 blocks = my_nb_data$parameters$C,
                                 sparsity = 0.1)
my_NB_sparse_high <- normal_block(data = my_data,
                                  blocks = my_nb_data$parameters$C,
                                  sparsity = 10)

## ----simulation1-NB-low-sparsity-plot-----------------------------------------
my_NB_sparse_low$plot_network()

## ----simulation1-NB-high-sparsity-plot----------------------------------------
my_NB_sparse_high$plot_network()

## ----simulation1-NB-changing-sparsity-----------------------------------------
my_NB_sparse <- normal_block(data = my_data,
                             blocks = my_nb_data$parameters$C,
                             sparsity = TRUE)

## ----simulation1-NB-changing-sparsity-plot------------------------------------
plot(my_NB_sparse)

## ----simulation1-NB-changing-sparsity-model-selection-------------------------
myNB_sparse_0.1   <- my_NB_sparse$get_model(0.1)
myNB_sparse_BIC   <- my_NB_sparse$get_best_model("BIC")

## ----simulation1-NB-sparse_unknown--------------------------------------------
my_NB_sparse_unknown <-  normal_block(data = my_data,
                                      blocks = 2:6,
                                      sparsity = TRUE)

## ----simulation1-NB-sparse_unknown-who-am-I-----------------------------------
my_NB_sparse_unknown$who_am_I

## ----simulation1-NB-sparse_unknown-plot---------------------------------------
plot(my_NB_sparse_unknown, "BIC")

## ----simulation1-NB-sparse_unknown-selection----------------------------------
my_NB_sparse_3     <- my_NB_sparse_unknown$get_model(3)
my_NB_sparse_3_0.1 <- my_NB_sparse_unknown$get_model(3, 0.1)

## ----simulation2-NB-----------------------------------------------------------
my_NB_zi <- normal_block(data = my_data_zi,
                         blocks = 4,
                         zero_inflation = TRUE)

## ----simulation2-NB-plot------------------------------------------------------
plot(my_NB_zi)

## ----simulation2-ARI----------------------------------------------------------
aricode::ARI(my_NB_zi$clustering,
             apply(my_nb_data_zi$parameters$C, 1, which.max))

