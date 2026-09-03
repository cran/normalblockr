## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----set-up-------------------------------------------------------------------
library(normalblockr)

## ----data-load----------------------------------------------------------------
data(brca_rppa)
dim(brca_rppa$expr)
table(brca_rppa$covariates$PAM50_SUBTYPE)

## ----NormalBlockData----------------------------------------------------------
Y            <- as.matrix(brca_rppa$expr)
X_subtype    <- model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates)
data_subtype <- NormalBlockData$new(Y, X_subtype)

## ----hclust-based-group, fig.width=7, fig.height=5----------------------------
hc_expr <- brca_rppa$expr |> scale() |> t() |> dist() |> hclust("ward.D2")
plot(hc_expr, labels = FALSE, hang = -1,
     main = "Hierarchical clustering of proteins (Ward, on scaled expression)",
     xlab = "proteins", sub = "")
rect.hclust(hc_expr, k = 6, border = "red")
group <- cutree(hc_expr, 6) |> normalblockr:::as_indicator()

## ----running-normal-block-known-group, fig.width=7, fig.height=5--------------
NB_prot_group <- normal_block(data_subtype, blocks = group)
plot(NB_prot_group)

## ----print-known-group--------------------------------------------------------
print(NB_prot_group)

## ----running-normal-block, fig.width=7, fig.height=5--------------------------
NB_prot_subtype <- normal_block(data_subtype, blocks = 1:30)

## ----plotting-criteria, fig.width=7, fig.height=5-----------------------------
NB_prot_subtype$plot(c("deviance", "BIC", "ICL"))

## ----model-selection----------------------------------------------------------
selected_NB <- NB_prot_subtype$get_best_model("ICL")
paste0("ICL selects ", selected_NB$q, " clusters.")

## ----refine-------------------------------------------------------------------
NB_prot_subtype$refine()

## ----plotting-criteria-refined, fig.width=7, fig.height=5---------------------
NB_prot_subtype$plot(c("deviance", "BIC", "ICL"))

## ----model-selection-refined--------------------------------------------------
selected_NB_refined <- NB_prot_subtype$get_best_model("ICL")
paste0("After refine(), ICL selects ", selected_NB_refined$q, " clusters.")

## ----sparsify-----------------------------------------------------------------
group_selected <- selected_NB_refined$clustering |> normalblockr:::as_indicator()
NB_prot_sparse <- normal_block(data_subtype, blocks = group_selected, sparsity = TRUE, control = NB_control(min_ratio=0.001))

## ----plotting-criteria-sparse, fig.width=7, fig.height=5----------------------
plot(NB_prot_sparse, c("EBIC", "deviance"))

## ----sparse-model-selection---------------------------------------------------
sparse_best <- NB_prot_sparse$get_best_model("EBIC")
paste0("BIC selects a penalty of ", round(sparse_best$sparsity, 4), ".")

## ----plot-network-sparse, fig.width=6, fig.height=6---------------------------
sparse_best$plot_network(output = "corrplot")

