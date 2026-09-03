## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----set-up-------------------------------------------------------------------
library(normalblockr)

## ----load-data----------------------------------------------------------------
data(onema)
dim(onema$biomass)

## ----zero-proportion----------------------------------------------------------
mean(onema$biomass == 0) # overall
range(colMeans(onema$biomass == 0)) # per-species range

## ----histogram-Y, fig.width = 6, fig.height = 4-------------------------------
Y_log <- log(1 + onema$biomass)
h <- hist(Y_log, breaks = 50, plot = FALSE)
ymax <- max(h$counts[-1]) * 1.15 # the zero bar (h$counts[1]) dwarfs every other one

plot(h, ylim = c(0, ymax), main = "", xlab = "log(1 + biomass)",
     col = "grey80", border = "white")
arrows(x0 = diff(h$breaks)[1] / 2, y0 = ymax * 0.93, y1 = ymax,
       length = 0.08, angle = 25, col = "grey30", lwd = 2)
text(x = diff(h$breaks)[1], y = ymax * 0.86,
     labels = paste0(round(mean(Y_log == 0) * 100), "% zeros, bar cut off"),
     pos = 4, col = "grey30")

## ----prepare-data-------------------------------------------------------------
X <- model.matrix(~ 1 + temperature_med, data = onema$covariates)
Y <- log(1 + onema$biomass)
data <- NormalBlockData$new(Y, X)

## ----fit-zi-------------------------------------------------------------------
out <- normal_block(data, blocks = 2:8, zero_inflation = TRUE)

## ----plot-zi------------------------------------------------------------------
out$plot(criteria = c("BIC", "EBIC", "deviance"))

## ----refine-zi----------------------------------------------------------------
out$refine()
out$plot(criteria = c("BIC", "EBIC", "deviance"))

## ----select-zi----------------------------------------------------------------
myModel <- out$get_best_model("EBIC")
myModel$q

## ----clusters-zi--------------------------------------------------------------
myModel$elements_per_cluster

## ----kappa-zi-----------------------------------------------------------------
kappa_hat <- myModel$model_par$kappa[1, ]
range(kappa_hat)

## ----kappa-check--------------------------------------------------------------
all.equal(kappa_hat, colMeans(onema$biomass == 0), check.attributes = FALSE)

## ----kappa-extremes-----------------------------------------------------------
sort(round(kappa_hat, 2), decreasing = TRUE)[1:5] # almost always structurally absent
sort(round(kappa_hat, 2))[1:5]                    # almost always present in some quantity

## ----sparsify-zi--------------------------------------------------------------
out_sp <- normal_block(data, blocks = myModel$q, sparsity = TRUE, zero_inflation = TRUE)
out_sp$plot(c("BIC", "EBIC"))
out_sp$plot(c("ICL", "deviance"))

## ----network-zi---------------------------------------------------------------
out_sp$get_best_model("ICL")$plot_network()

