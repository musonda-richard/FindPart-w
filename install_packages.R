packages <- c("tidyverse","lubridate","igraph", "ggraph", "ggrepel", "grid")
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
lapply(packages, library, character.only = TRUE)
