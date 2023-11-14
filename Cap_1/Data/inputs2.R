
#header #################################################################################
#'inputs2.R'

#Title: Support file
#Project ID: pid
#Client: client
#Author: <Eduardo> <Costa>, 

#Description: Support file with inputs for the univariate sensitivity analysis.

#Start date: date
#Last Update: {6:date}

#R version: r.version
#Scriptversion: version

#Dependencies
#<-Downstream
# + None

#->Upstream
# + model_carcas_univariate.R

#Input:
# +None

#Output:
# +None

#Peer reviewer(s)

#Please ensure directories are relative. Please comment code sufficiently.

#Script start#############################################################################


x<-c(-1,-0.85,-0.75,-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,0.4,0.5,0.75,0.85,1)

u<-function(x, base,min,max){
  
  ifelse(x>=0,base-(base-max)*x,
       base+(base-min)*x)
  }


## Inputs

#Transf prob from environm to carcass/organ
ak_c=0.0017 #knife-carcass

#hand-carcass
ah_cbase=0.0021
ah_cmin=0
ah_cmax=1

ah_c=list(rep(ah_cbase,length(x)),
          rep(ah_cbase,length(x)),
          rep(ah_cbase,length(x)),
          rep(ah_cbase,length(x)),
          rep(ah_cbase,length(x)),
          rep(ah_cbase,length(x)),
          u(x,base = ah_cbase,min=ah_cmin,max=ah_cmax)
          
          )

  
ag_c=0.0017 #hook-carcass
ak_o=0.0017 #knife-organ


#hand-organ
ah_obase=0.0021
ah_omin=0
ah_omax=1

ah_o=list(rep(ah_obase,length(x)),
          rep(ah_obase,length(x)),
          rep(ah_obase,length(x)),
          rep(ah_obase,length(x)),
          rep(ah_obase,length(x)),
          u(x,base = ah_obase,min = ah_omin,max=ah_omax),
          rep(ah_obase,length(x))
          
          )
  


ag_o=0.0017 #hook-organ


#Transf prob from carcass/organ to environment
bc_k=0.0017  #carcass-knife
bc_h=0.031   #carcass-hand
bc_g=0.0017   #carcass-hook
bo_k=0.0017  #oragan-knife


#organ-hand
bo_hbase=0.031
bo_hmin=0
bo_hmax=1

bo_h=list(rep(bo_hbase,length(x)),
          rep(bo_hbase,length(x)),
          rep(bo_hbase,length(x)),
          rep(bo_hbase,length(x)),
          u(x,base =bo_hbase,min = bo_hmin,max=bo_hmax),
          rep(bo_hbase,length(x)),
          rep(bo_hbase,length(x))
        )
  
  
bo_g=0.0017 #organ-hook


#Prob. of destruiction in cacass/organ
dc=0 #destruct carcass
do=0 #destruct organ

##Touch area
ae_k=10  #knife
ae_h=150 #hand
ae_g=1   #hook

##Number of contacts
co=seq(from=0, to=22,by=1)

##Prob of touches
p_h=c(0.0089,0.0089,0.01785,0.02678,0.1517,0.116,0.1428,0.1339,0.1607,0.0464,
      0.0982,0.0089,0.01785,0,0.01785,0.008929,0.01785,0,0,0.008929,0.008929,0,0)

##Prob of hook
p_g=c(0,0.0166,0.4833,0.433,0.033,0.0166,0.0166,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

##Prob of cuts
p_k=c(0,0.008,0.051,0.034,0.144,0.161,0.169,0.161,0.059,0.068,0.051,0,0.034,0,0.017,0.025,0,0.008,0,0.008,0,0,0)

## Prob of change knife
ck=0.9

# initial contamination log10 cfu/cm^2

mu_base=-5.4
mu_min=-15
mu_max=2

mu=list(u(x,base=mu_base,min = mu_min,max=mu_max),
        rep(mu_base,length(x)),
        rep(mu_base,length(x)),
        rep(mu_base,length(x)),
        rep(mu_base,length(x)),
        rep(mu_base,length(x)),    
        rep(mu_base,length(x))
        )

sig_base=2.2
sig_min=0
sig_max=4


sig=list(rep(sig_base,length(x)),
         u(x,base=sig_base,min = sig_min,max=sig_max),
         rep(sig_base,length(x)),
         rep(sig_base,length(x)),
         rep(sig_base,length(x)),
         rep(sig_base,length(x)),
         rep(sig_base,length(x))
  
)


#initial prevalence carcass
prev=1

#initial prevalence LN
alpha_obase=34
alpha_omin=1
alpha_omax=140

alpha_o=list(rep(alpha_obase,length(x)),
             rep(alpha_obase,length(x)),
             u(x,base=alpha_obase,min = alpha_omin,max=alpha_omax),
             rep(alpha_obase,length(x)),
             rep(alpha_obase,length(x)),
             rep(alpha_obase,length(x)),
             rep(alpha_obase,length(x))
)


beta_o<-107

#concentration in organ cfu/cm^2
min=0.01
mod=1

max_base=100
max_min=10
max_max=1000
max=list(rep(max_base,length(x)),
         rep(max_base,length(x)),
         rep(max_base,length(x)),
         u(x,base=max_base,min = max_min,max=max_max),
         rep(max_base,length(x)),
         rep(max_base,length(x)),
         rep(max_base,length(x))
  
)
  

##Crating a function to recognize integers
check.integer <- function(x) {
  x == round(x)
}


##Creating objects for simulation

K=numeric()
H=numeric()
G=numeric()

j_k=numeric()
j_h=numeric()
j_g=numeric()

jo_k=numeric()
jo_h=numeric()
jo_g=numeric()

jc_k=numeric()
jc_h=numeric()
jc_g=numeric()

C_S1=numeric()
C_kS=numeric()
C_hS=numeric()
C_gS=numeric()

C_k1=numeric()
C_h1=numeric()
C_g1=numeric()

C_k=numeric()
C_h=numeric()
C_g=numeric()


C_k2=numeric()
C_h2=numeric()
C_g2=numeric()

C_S=numeric()
C=numeric()

O_kS=numeric()
O_hS=numeric()
O_gS=numeric()

O_k1=numeric()
O_h1=numeric()
O_g1=numeric()


O_k=numeric()
O_h=numeric()
O_g=numeric()

O_k2=numeric()
O_h2=numeric()
O_g2=numeric()

O_S=numeric()
O=numeric()

C_S_conc=numeric()
C_conc=numeric()

O_S_conc=numeric()

O_conc=numeric()

prev_CS=numeric()

prev_C=numeric()

prev_OS=numeric()

prev_O=numeric()


logCS=numeric()
logC=numeric()

logOS=numeric()
logO=numeric()

O_S1=numeric()  

O_kS=numeric()
O_hS=numeric()
O_gS=numeric()

O_k1=numeric()
O_h1=numeric()
O_g1=numeric()


O_k=numeric()
O_h=numeric()
O_g=numeric()


O_k2=numeric()
O_h2=numeric()
O_g2=numeric()

O_S=numeric()
O=numeric()

C_S_conc=numeric()
C_conc=numeric()


O_S_conc=numeric()
O_conc=numeric()

prev_CS=numeric()

prev_C=numeric()

prev_OS=numeric()

prev_O=numeric()


logCS=numeric()
logC=numeric()

logOS=numeric()
logO=numeric()


me_CS=numeric()
me_C=numeric()
me_OS=numeric()
me_O=numeric()

sd_CS=numeric()
sd_C=numeric()
sd_OS=numeric()
sd_O=numeric()

prevalence_CS=numeric()
prevalence_C=numeric()
prevalence_OS=numeric()
prevalence_O=numeric()


mu_mu_c=vector(mode = "list", length = 7)
mu_p_c=vector(mode = "list", length = 7)



names<-c("Mean[c_i(S-1)]",
         "SD[c_i(S-1)]",
         "Prev[O_i(S-1)]",
         "Max[O_i(S-1)]",
         "b_OH",
         "a_HO",
         "a_HC")
