# EcoState

[![](https://www.r-pkg.org/badges/version/ecostate)](https://cran.r-project.org/package=ecostate)
[![](https://cranlogs.r-pkg.org/badges/ecostate)](https://cran.r-project.org/package=ecostate)
[![](https://cranlogs.r-pkg.org/badges/grand-total/ecostate)](https://cran.r-project.org/package=ecostate)
[![Documentation](https://img.shields.io/badge/documentation-ecostate-orange.svg?colorB=E91E63)](https://james-thorson-noaa.github.io/ecostate/)
[![Codecov test coverage](https://codecov.io/gh/James-Thorson-NOAA/ecostate/branch/add_dsem/graph/badge.svg)](https://app.codecov.io/gh/James-Thorson-NOAA/ecostate/tree/add_dsem)

Package _ecostate_ fits a state-space mass-balance model _EcoState_ intended for aquatic ecosystems, using mass-balance equations matching those from Ecopath and dynamical equations matching Ecosim.  Unlike the Ecopath-with-Ecosim (_EwE_) package, _ecostate_ fits both biological parameters (e.g., equilibrium biomass and predator-prey vulnerability) and measurement parameters (e.g., catchability coefficients) via fit to time-series data.  _ecostate_ also estimates additional process errors representing nonstationarity in growth efficiency, ecotrophic efficient, migration, or other unmodeled processes.  These process errors allow biomass patterns to closely match available data, so that resulting consumption (and associated productivity and mortality rates) can accurately be conditioned upon any residual patterns.     

## Installation

_ecostate_ can be installed from GitHub using:

``` r
library(remotes)
install_github( "James-Thorson-NOAA/ecostate" )
```

Or to access vignettes from your R session, please instead use:

``` r
remotes::install_github( "James-Thorson-NOAA/ecostate",
                          build_vignettes = TRUE )
browseVignettes("ecostate")
```

# More details 

For more background, please read:

### For state-space mass-balance equilibrium and biomass-dynamic modelling
Thorson, J.  Kristensen, K., Aydin, K., Gaichas, S., Kimmel, D.G., McHuron, E.A., Nielsen, J.N., Townsend, H., Whitehouse, G.A. (2025) The Benefits of Hierarchical Ecosystem Models: Demonstration Using EcoState, a New State-Space Mass-Balance Mode.  Fish and Fisheries, 26(2): 203-218. https://doi.org/10.1111/faf.12874

### For age-structured dynamics and bottom-up linkages to animal growth:
Thorson, J. T., Aydin, K. H., Cheng, M., Dias, B. S., Kimmel, D. G., & Kristensen, K. (2025). Bottom-up interactions in age-structured stock assessment and state-space mass-balance modelling. http://doi.org/10.1111/faf.70016

### For original development of ecopath for mass-balance equilibria:
Polovina, J. J. (1984). Model of a coral reef ecosystem. Coral Reefs, 3(1), 1–11. https://doi.org/10.1007/BF00306135

### For development of ecosim for deterministic biomass and age-structured dynamics:
Christensen, V., & Walters, C. J. (2004). Ecopath with Ecosim: Methods, capabilities and limitations. Ecological Modelling, 172(2), 109–139. https://doi.org/10.1016/j.ecolmodel.2003.09.003

Walters, C., Pauly, D., Christensen, V., & Kitchell, J. F. (2000). Representing Density Dependent Consequences of Life History Strategies in Aquatic Ecosystems: EcoSim II. Ecosystems, 3(1), 70–83. https://doi.org/10.1007/s100210000011

# NOAA Enterprise GitHub disclaimer
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

