
## Inputs
#Transf prob from carcass/organ to environment
bc_k=c(0.0017,0,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017)  #carcass-knife
bc_h=c(0.0310,0,0.0310,0.0310,0.0310,0.1,0.1,0.1,0.1)   #carcass-hand
bc_g=c(0.0017,0,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017)  #carcass-hook
bo_k=c(0.0017,0,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017)  #oragan-knife
bo_h=c(0.0310,0,0.0310,0.0310,0.0310,0.0310,0.0310,0.0310,0.0310)  #organ-hand
bo_g=c(0.0017,0,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017)  #organ-hook

#Transf prob from environm to carcass/organ
ak_c=c(0.0017,0,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017) #knife-carcass
ah_c=c(0.0021,0,0.0021,0.0021,0.0021,0.1,0.1,0.1,0.1) #hand-carcass
ag_c=c(0.0017,0,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017) #hook-carcass
ak_o=c(0.0017,0,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017) #knife-organ
ah_o=c(0.0021,0,0.0021,0.0021,0.0021,0.0021,0.0021,0.0021,0.0021) #hand-organ
ag_o=c(0.0017,0,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017,0.0017) #hook-organ

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
mu=c(-5.4,-5.4,-5.4,-5.4,0,-3,-5.4,-5.4,-5.4)
sig=c(2.2,2.2,2.2,2.2,0,2.2,3.3,2.2,2.2)

#initial prevalence carcass
prev=c(1,1,0,1,0,1,1,1,1)

#initial prevalence LN
alpha_o<-c(34,34,34,34,34,34,34,34,73)
beta_o<-c(107,107,107,107,107,107,107,107,1)

#concentration in organ cfu/cm^2
min=c(0.01,0.01,0,0,0.01,0.01,0.01,0.1,0.01)
mod=c(1,1,0,0,1,1,1,10,1)
max=c(100,100,0,0,100,100,100,1000,100)


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

me_CS1=numeric()
me_C1=numeric()

sd_CS=numeric()
sd_C=numeric()
sd_OS=numeric()
sd_O=numeric()

prevalence_CS=numeric()
prevalence_C=numeric()
prevalence_OS=numeric()
prevalence_O=numeric()

suma=list()


names<-c("Baseline",
         "Null transfer",
         "Null contaminaiton",
         "Only Carcass (S-1)",
         "Only LN (S-1)",
         "Hand influence High Mean Carcass (S-1)",
         "Hand influence High SD Carcass (S-1)",
         "Hand Influence high LN CFU (S-1)",
         "Hand Influence high prev LN (S-1)")



