[![DOI](https://zenodo.org/badge/DOI/10.1111/risa.12753.svg)](https://doi.org/10.1111/risa.12753)

:pig: :hand: :knife: :mag: :cut_of_meat:


# Title: 

A Stochastic Model to Assess the Effect of Meat Inspection Practices on the Contamination of the Pig Carcasses



# Authors: 

Eduardo de Freitas Costa, Luis Gustavo Corbellini, Ana Paula Serafini Poeta Silva and Maarten Nauta



## The paper is available online:

https://onlinelibrary.wiley.com/doi/abs/10.1111/risa.12753

Costa, de F. C., Corbellini, L. G., Poeta Silva, A. P. S., Nauta, M. (2016). A Stochastic Model to Assess the Effect of Meat Inspection Practices on the Contamination of the Pig Carcasses. **Risk Analysis**, V. 37 (10), p, 1849-1864. DOI: 10.1111/risa.12753.

## This is an @Risk file for the cros-contamination of *Salmonella* spp. during inspeciton of carcasses of pigs in Brazil



# Instrucitons to reproduce the results from the paper "A Stochastic Model to Assess the Effect of Meat Inspection Practices on the Contamination of the Pig Carcasses".

## Using R or RStudio.

### 1 Download the files to a local or cloud folder. The best option is to download all files zipped;

### 2 Unzip and open the R project "Doutorado.Rproj" (it is highly recommended to use RStudio);

### 3 In R, you can open the file\ 

  + "main.R" \ 
  
### 4 Load the packages before running the model;

### 5 After that, choose if the baseline, univariate or multivariate models will be run. All functions follow the same:  *model(ncarc=,nrepl=)*, where *ncarc* is the number of carcasses considered, and *nrepl* is the number of stochastic iterations;

### Note that two folders were created into the folder *Cap_1*:
  - **Output**: Stores the results of the baseline, univariable, and multivariate scenarios simulation;
  
  - **Figures**: Stores the graphs for the *Salmonella* sp. concentration on carcasses surfaces before and after inspection for all scenarios, including the baseline;
  
## Using @Risk

### 1 The user needs the add in @Risk installed. After that, copy and open the file\
 + "inspection_model_12.xls" (available in the tab "Files" into the *Cap_1\Scripts* folder)\ 
 
### 2 Execute the stochastic simulation according to the @Risk command. The results will be available in the outputs tabs. 

