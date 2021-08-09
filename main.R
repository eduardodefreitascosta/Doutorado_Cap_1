

# Installing and loading packages/libraries


packages<-c("here","triangle","ggplot2","psych")

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


#Run the script

source(here("Cap_1","Scripts","model_carcas.R"))