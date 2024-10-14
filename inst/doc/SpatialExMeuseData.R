## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(kableExtra)

## ----echo = FALSE, results='hide', message=FALSE, warning=FALSE---------------
library(GGally)
# refactor from sp to sf when meuse dataset available through sf
library(ggspatial)
library(spData)
library(sf)
library(cowplot)
library(giscoR)
library(wesanderson)
library(magrittr)
library(kableExtra)
library(knitr)
library(ggplot2)
library(inlabru)
library(fmesher)

## -----------------------------------------------------------------------------
# load data
library(sp)
data("meuse")
data("meuse.riv")
data("meuse.grid")

## ----sample, fig.height = 4, fig.width = 6, fig.cap = "Sample locations of heavy metal concentration levels (ppm) of four metals in the topsoil of the Meuse river floodplain, Netherlands.",  echo = FALSE----
# map study area
# read in and project spatial data
m.river <- as.data.frame(meuse.riv[meuse.riv[, 2] > 329000 & meuse.riv[, 2] < 334200, ])
colnames(m.river) <- c("x", "y")
m.river <- sf::st_as_sf(m.river,
  coords = c("x", "y"),
  crs = 28992
) %>%
  dplyr::summarise(geometry = sf::st_combine(geometry)) %>%
  sf::st_cast("POLYGON")

m.grid <- sf::st_as_sf(meuse.grid,
  coords = c("x", "y"),
  crs = 28992
)
m.data <- sf::st_as_sf(meuse,
  coords = c("x", "y"),
  crs = 28992
)

mRiv <- sf::st_transform(m.river, crs = 4326)
mGrid <- sf::st_transform(m.grid, crs = 4326)
mDat <- sf::st_transform(m.data, crs = 4326)

# create map breaks and labels
xbreaks <- seq(5.72, 5.76, .01)
xlabels <- unlist(sapply(xbreaks, function(x) paste(x, "°E")))
ybreaks <- seq(50.95, 51.00, .01)
ylabels <- unlist(sapply(ybreaks, function(x) paste(x, "°N")))
# create reference map
Europe <- giscoR::gisco_get_countries(region = "Europe")
r.bb <- sf::st_bbox(mRiv)
r.bb[1] <- r.bb[1] - 0.5
r.bb[2] <- r.bb[2] - 0.3
r.bb[3] <- r.bb[3] + 0.5
r.bb[4] <- r.bb[4] + 0.3
rbox <- sf::st_as_sfc(r.bb)
ggm1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = Europe, fill = "cornsilk") +
  ggplot2::geom_sf(data = rbox, fill = NA, color = "red", size = 1.2) +
  coord_sf(
    xlim = c(-5, 20),
    ylim = c(46.5, 56.5)
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "#c9f4fe"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1.5
    )
  )

# create site map
ggm2 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = mDat) +
  theme_classic() +
  geom_sf(data = mRiv, fill = "lightblue") +
  ggspatial::annotation_north_arrow(
    which_north = "grid", height = unit(1, "cm"), width = unit(1, "cm"),
    pad_y = unit(4, "cm")
  )
cowplot::ggdraw() + cowplot::draw_plot(ggm2) + cowplot::draw_plot(ggm1, x = .18, y = 0.7, width = 0.35, height = 0.23)

## ----pairs-expl, fig.height = 4, fig.width = 6,  fig.cap = "Pairs plots of the concentration levels (ppm) of four heavy metals measured in the topsoil of the Meuse river floodplain, Netherlands.", echo = FALSE----
# Visualize heavy metal data
ggpairs(meuse[, 3:6])

## ----warning = FALSE----------------------------------------------------------
library(clustTMB)
mod.gauss <- clustTMB(response = meuse[, 3:6], G = 3, covariance.structure = "VVV")

## -----------------------------------------------------------------------------
mod.ln <- clustTMB(
  response = meuse[, 3:6],
  family = lognormal(link = "identity"),
  G = 3, covariance.structure = "VVV"
)

## -----------------------------------------------------------------------------
library(fmesher)
loc <- meuse[, 1:2]
Bnd <- fmesher::fm_nonconvex_hull(as.matrix(loc), convex = 200)
meuse.mesh <- fmesher::fm_mesh_2d(as.matrix(loc),
  max.edge = c(300, 1000),
  boundary = Bnd
)

## ----fig1, fig.height = 3, fig.width = 4, warning = FALSE---------------------
ggplot() +
  inlabru::gg(meuse.mesh) +
  geom_point(mapping = aes(x = loc[, 1], y = loc[, 2], size = 0.5), size = 0.5) +
  theme_classic()

## ----warning = FALSE----------------------------------------------------------
# convert coordinates to a spatial point data frame
Loc <- sf::st_as_sf(loc, coords = c("x", "y"))

# define spatial prediction coordinates
data("meuse.grid")
Meuse.Grid <- sf::st_as_sf(meuse.grid, coords = c("x", "y"))
mod.ln.sp <- clustTMB(
  response = meuse[, 3:6],
  family = lognormal(link = "identity"),
  gatingformula = ~ gmrf(0 + 1 | loc),
  G = 4, covariance.structure = "VVV",
  spatial.list = list(loc = Loc, mesh = meuse.mesh),
  projection.dat = Meuse.Grid
)

## -----------------------------------------------------------------------------
# Estimated fixed parameters
summary(mod.ln.sp$sdr, "fixed")
# Minimum negative log likelihood
mod.ln.sp$opt$objective

## ----BIC, echo = FALSE--------------------------------------------------------
df <- data.frame(
  family = c(rep("lognormal", 2), rep("Gaussian", 2)),
  space = c(1, 0, 1, 0),
  clusters = c(4, 4, 6, 3),
  BIC = c(4805, 4861, 4861, 4910)
)
kable_styling(
  kable(df,
    booktabs = TRUE,
    caption = "Optimum cluster size and BIC scores for lognormal and Gaussian models with (1) and without (0) spatial random effects in the gating network."
  ),
  full_width = TRUE
)

## ----pairs, fig.height = 4, fig.width = 6, fig.cap = "Pairs plots of the concentration levels (ppm) of four heavy metals measured in the topsoil of the Meuse river floodplain, Netherlands. Colors represent the four clusters estimated in the spatial lognormal FMM.", echo = FALSE----
# Visualize heavy metal data
Cluster <- factor(mod.ln.sp$report$classification)
ggpairs(meuse[, 3:6], aes(color = Cluster, alpha = 90.5)) +
  scale_color_manual(values = wes_palette("Moonrise2", 4)) +
  scale_fill_manual(values = wes_palette("Moonrise2", 4))

## ----pred, fig.height = 4, fig.width = 6, fig.cap = "Predicted cluster distribution of heavy metal concentration levels (ppm) of four metals in the topsoil of the Meuse river floodplain, Netherlands. Colors represent the four clusters estimated in the spatial lognormal FMM.", echo = FALSE, warning = FALSE----
Cluster <- factor(mod.ln.sp$report$Class_pred)

# create site map
ggm2 <- ggplot2::ggplot(mGrid) +
  ggplot2::stat_sf_coordinates(aes(color = Cluster), geom = "point") +
  theme_classic() +
  scale_color_manual(values = wes_palette("Moonrise2", 4)) +
  geom_sf(data = mRiv, fill = "lightblue") +
  ggspatial::annotation_north_arrow(
    which_north = "grid", height = unit(1, "cm"), width = unit(1, "cm"),
    pad_y = unit(4, "cm")
  )
cowplot::ggdraw() + cowplot::draw_plot(ggm2) + cowplot::draw_plot(ggm1, x = .18, y = 0.7, width = 0.35, height = 0.23)

