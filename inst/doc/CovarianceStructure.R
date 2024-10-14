## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(kableExtra)

## ----table, echo = FALSE, warnings = FALSE------------------------------------
tbl <- data.frame(
  "Covariance" = c("Spatial GMRF", "AR(1)", "Rank Reduction", "Spatial Rank Reduction"),
  "Notation" = c("gmrf", "ar1", "rr(random = H)", "rr(spatial = H)"),
  "No. of Parameters" = c("2", "2", "JH - (H(H-1))/2", "1 + JH - (H(H-1))/2"),
  "Data requirements" = c("spatial coordinates", "unit spaced levels", "", "spatial coordinates")
)
kbl(tbl, booktabs = TRUE)

## ----spatial example, include = FALSE-----------------------------------------
library(clustTMB)
# refactor from sp to sf when meuse dataset available through sf
library(sp) # currently require sp to load meuse dataset
data("meuse")
library(fmesher)

## ----meuse mesh---------------------------------------------------------------
loc <- meuse[, 1:2]
Bnd <- fmesher::fm_nonconvex_hull(as.matrix(loc), convex = 200)
meuse.mesh <- fmesher::fm_mesh_2d(as.matrix(loc),
  max.edge = c(300, 1000),
  boundary = Bnd
)

## ----fig1, fig.height = 3, fig.width = 5, echo = FALSE------------------------
library(ggplot2)
library(inlabru)
ggplot() +
  gg(meuse.mesh) +
  geom_point(mapping = aes(x = loc[, 1], y = loc[, 2], size = 0.5), size = 0.5) +
  theme_classic()

## ----set up model-------------------------------------------------------------
Loc <- sf::st_as_sf(loc, coords = c("x", "y"))
mod <- clustTMB(
  response = meuse[, 3:6],
  family = lognormal(link = "identity"),
  gatingformula = ~ gmrf(0 + 1 | loc),
  G = 4, covariance.structure = "VVV",
  spatial.list = list(loc = Loc, mesh = meuse.mesh)
)

## ----view results-------------------------------------------------------------
# Estimated fixed parameters
mod$opt$par
# Minimum negative log likelihood
mod$opt$objective

