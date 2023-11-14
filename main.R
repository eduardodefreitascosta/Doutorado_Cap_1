

# Installing and loading packages/libraries
packages<-c("here","triangle","ggplot2","psych","scales","gridExtra","tidyverse")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))  


# Creating directories
dir.create(here("Cap_1","Output"),showWarnings = F)
dir.create(here("Cap_1","Figures"),showWarnings = F)


#Run the baseline model
source(here("Cap_1","Scripts","model_carcas.R"))
model(ncarc=500, nrepl=5000)


#Run the univariable sensitivity model
source(here("Cap_1","Scripts","model_carcas_univariable.R"))
model_uni(ncarc=500, nrepl=5000)


#Run the multivariable sensitivity model
source(here("Cap_1","Scripts","model_carcas_multi.R"))
model(ncarc=50, nrepl=50)


