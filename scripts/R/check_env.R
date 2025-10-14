dir.create(
  path = Sys.getenv("R_LIBS_USER"),
  showWarnings = FALSE,
  recursive = TRUE
)

packages <- c(
  "tidyverse", # nolint
  "patchwork",
  "arsenal",
  "finalfit",
  "ggsci",
  "reticulate",
  "ggalluvial",
  "ggrepel",
  "ggnewscale",
  "scico",
  "jsonlite",
  "base64enc"
)

# Install packages that are not already installed
install.packages(
  setdiff(packages, rownames(installed.packages())),
  lib = Sys.getenv("R_LIBS_USER"),
  repos = "https://cran.rstudio.com/"
)