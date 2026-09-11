## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----set-up-------------------------------------------------------------------
library(normalblockr)

## ----data-load----------------------------------------------------------------
data(brca_rppa)
Y <- as.matrix(brca_rppa$expr)
X <- model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates)
nb_data <- NormalBlockData$new(Y, X)
dim(Y)
table(brca_rppa$covariates$PAM50_SUBTYPE)

## ----known-clustering---------------------------------------------------------
go_term  <- factor(brca_rppa$gene_annotation$go_bp_term)
C_go     <- model.matrix(~ 0 + go_term)
NB_go    <- normal_block(nb_data, blocks = C_go, model = "mean",
                         control = NB_control(verbose = FALSE))
NB_go

## ----known-clustering-fitted, fig.width=6, fig.height=4-----------------------
dim(coef(NB_go))
plot(Y, fitted(NB_go), pch = ".", xlab = "observed", ylab = "fitted")
abline(0, 1, col = "red")

## ----collection---------------------------------------------------------------
NB_means <- normal_block(nb_data, blocks = seq(5, 120, by = 5), model = "mean",
                         control = NB_control(verbose = FALSE))

## ----collection-criteria, fig.width=7, fig.height=4---------------------------
NB_means$plot(c("deviance", "ICL"))

## ----collection-selection-----------------------------------------------------
selected <- NB_means$get_best_model("ICL")
paste0("ICL selects ", selected$q, " clusters.")

## ----clusters-----------------------------------------------------------------
table(selected$clustering)

## ----cluster-profiles, fig.width=7, fig.height=4------------------------------
profiles <- coef(selected)
rownames(profiles) <- levels(brca_rppa$covariates$PAM50_SUBTYPE)
colnames(profiles) <- paste0("cluster ", seq_len(ncol(profiles)))
image(seq_len(nrow(profiles)), seq_len(ncol(profiles)), profiles,
      axes = FALSE, xlab = "", ylab = "", col = hcl.colors(20, "RdBu", rev = TRUE))
axis(1, seq_len(nrow(profiles)), rownames(profiles), las = 2, cex.axis = .7)
axis(2, seq_len(ncol(profiles)), colnames(profiles), las = 2, cex.axis = .7)

## ----covariance-shapes--------------------------------------------------------
shapes <- c("diagonal", "spherical", "full")
fits   <- lapply(shapes, function(s)
  normal_block(nb_data, blocks = selected$q, model = "mean",
               control = NB_control(verbose = FALSE, noise_covariance = s)))
data.frame(
  covariance = shapes,
  nb_param   = sapply(fits, `[[`, "nb_param"),
  loglik     = round(sapply(fits, `[[`, "loglik"), 1),
  BIC        = round(sapply(fits, `[[`, "BIC"), 1)
)

## ----sparse-fit---------------------------------------------------------------
C_selected <- model.matrix(~ 0 + factor(selected$clustering))
NB_sparse  <- normal_block(nb_data, blocks = C_selected, sparsity = 0.4,
                           model = "mean", control = NB_control(verbose = FALSE))
NB_sparse$model_par$Omega |> dim()
paste0(NB_sparse$n_edges, " edges out of ", choose(ncol(Y), 2), " possible ones.")

## ----sparse-network, fig.width=6, fig.height=6--------------------------------
NB_sparse$plot_network(output = "corrplot")

