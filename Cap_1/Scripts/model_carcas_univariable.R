
#header #################################################################################
#'model_carcass_univariable.R'

#Title: Stochastic modeling for cross-contamination carcasses during inspection
#Project ID: pid
#Client: client
#Author: <Eduardo> <Costa>, 

#Description: Base line and univariable sensitivity analysis for cross-contamination during 
# pigs carcassess inspection in BR

#Start date: date
#Last Update: {6:date}

#R version: r.version
#Scriptversion: version

#Dependencies
#<-Downstream
# +inputs2.r

#->Upstream
# +main.r

#Input:
# +Data/inputs2.R

#Output:
# +Uni_scenarios/ .txt
# +Figures/Uni_scenarios/ .png

#Peer reviewer(s)

#Please ensure directories are relative. Please comment code sufficiently.

#Script start#############################################################################

model_uni<-function(ncarc, nrepl){
  
  
  #Calling inputs
  source(here("Cap_1","Data","inputs2.R"))
  
  ######################################
  ##Begin of the external list looping##
  ######################################

  for(k in 1:7){
    tempo_inicial1 = proc.time()
    cat("\n\n scenario: ", k, "\n")
    
    ######################################
    ##Begin of the internal list looping##
    ######################################
    
    for(w in 1:length(x)){
      tempo_inicial1 = proc.time()
      cat("\n\n scenario: ", w, "\n")
      
    
    ###################################
    ##Begin of the stochastic looping##
    ###################################
    
    for(j in 1:nrepl){
      tempo_inicial1 = proc.time()
      cat("\n\n nrepl: ", j, "\n")
      
      ####Initial states
      ###environment
      K[1]=0
      H[1]=0
      G[1]=0
      
      j_k[1]=round(sample(co,size=1,prob=p_k,replace=T),0) #total number of cuts
      j_h[1]=round(sample(co,size=1,prob=p_h,replace=T),0) #total number of handling
      j_g[1]=round(sample(co,size=1,prob=p_g,replace=T),0) #total number of hook
      
      
      jo_k[1]=rbinom(1,j_k[1],0.5) #number of cuts organs
      jo_h[1]=rbinom(1,j_h[1],0.5) #number of handling organs
      jo_g[1]=rbinom(1,j_g[1],0.5) #number of hooks organs
      
      
      jc_k[1]=j_k[1]-jo_k[1] #number of cuts carcass
      jc_h[1]=j_h[1]-jo_h[1] #number of handling carcass
      jc_g[1]=j_g[1]-jo_g[1] #number of hook carcass
      
      
      
      
      ###Carcass
      
      ##S-1
      C_S1[1]=rpois(1,10^rnorm(1,mu[[k]][w],sig[[k]][w]))*rbinom(1,1,prev) #carcass concentration surface
      
      C_kS[1]=(rpois(1,C_S1[1]*ae_k*jc_k[1])) #cfu in knife area
      C_hS[1]=(rpois(1,C_S1[1]*ae_h*jc_h[1])) #cfu in hand area
      C_gS[1]=(rpois(1,C_S1[1]*ae_g*jc_g[1])) #cfu in hook area
      
      ##S
      #Detection
      C_k1[1]=rbinom(1,C_kS[1],1-dc) #cfu after detection knife area
      C_h1[1]=rbinom(1,C_hS[1],1-dc) #cfu after detection hand area
      C_g1[1]=rbinom(1,C_gS[1],1-dc) #cfu after detection hook area
      
      #Loose to environment
      C_k[1]=rbinom(1,C_k1[1],(1-bc_k)^jc_k[1]) #cfu 1-loose to knife
      C_h[1]=rbinom(1,C_h1[1],(1-bc_h)^jc_h[1]) #cfu 1-loose to hand
      C_g[1]=rbinom(1,C_g1[1],(1-bc_g)^jc_g[1]) #cfu 1-loose to hook
      
      #Gain from environment
      C_k2[1]= K[1]-rbinom(1,K[1],(1-ak_c)^jc_k[1]) #gain form knife
      C_h2[1]= H[1]-rbinom(1,H[1],(1-ah_c[[k]][w])^jc_h[1]) #gain form hands 
      C_g2[1]= G[1]-rbinom(1,G[1],(1-ag_c)^jc_g[1]) #gain form hook
      
      C_S[1]=C_kS[1]+C_hS[1]+C_gS[1] #Surface before inspection
      C[1]=C_k[1]+C_k2[1]+C_h[1]+C_h2[1]+C_g[1]+C_g2[1] #Surface after insepction
      
      
      ###Limph nodes
      ##S-1
      O_S1[1]=rpois(1,rtriangle(n=1,a=min,b=max[[k]][w],c=mod))*rbinom(1,1,alpha_o[[k]][w]/(alpha_o[[k]][w]+beta_o)) # LN concentration
      
      
      O_kS[1]=rpois(1,O_S1[1]*ae_k*jo_k[1]) #cfu in knife area
      O_hS[1]=rpois(1,O_S1[1]*ae_h*jo_h[1]) #cfu in hand area
      O_gS[1]=rpois(1,O_S1[1]*ae_g*jo_g[1]) #cfu in hook area
      
      ##S
      #Detection
      O_k1[1]=rbinom(1,O_kS[1],1-do)  #cfu after detection knife area
      O_h1[1]=rbinom(1,O_hS[1],1-do)  #cfu after detection hand area
      O_g1[1]=rbinom(1,O_gS[1],1-do)  #cfu after detection hook area
      
      #Loose environment
      O_k[1]=rbinom(1,O_k1[1],(1-bo_k)^jo_k[1])  #cfu 1-loose to knife
      O_h[1]=rbinom(1,O_h1[1],(1-bo_h[[k]][w])^jo_h[1])   #cfu 1-loose to hand
      O_g[1]=rbinom(1,O_g1[1],(1-bo_g)^jo_g[1])   #cfu 1-loose to hook
      
      #Gain from environment
      O_k2[1]= (K[1]-C_k2[1])-rbinom(1,(K[1]-C_k2[1]),(1-ak_o)^jo_k[1]) #gain form knife
      O_h2[1]= (H[1]-C_h2[1])-rbinom(1,(H[1]-C_h2[1]),(1-ah_o[[k]][w])^jo_h[1]) #gain form hands 
      O_g2[1]= (G[1]-C_g2[1])-rbinom(1,(G[1]-C_g2[1]),(1-ag_o)^jo_g[1]) #gain form hook
      
      O_S[1]=O_kS[1]+O_hS[1]+O_gS[1] #LN before inspection
      O[1]=O_k[1]+O_k2[1]+O_h[1]+O_h2[1]+O_g[1]+O_g2[1] #LN after inspection
      
      
      
      ##cfu/cm^2
      #Crcass
      C_S_conc[1]=ifelse((jc_k[1]+jc_h[1]+jc_g[1])==0,C_S1[1],  #concentration before inspection
                         C_S[1]/((jc_k[1]*ae_k)+(jc_h[1]*ae_h)+(jc_g[1]*ae_g))
      )
      
      C_conc[1]=ifelse((jc_k[1]+jc_h[1]+jc_g[1])==0,C_S1[1],  #concentration after inspection
                       C[1]/((jc_k[1]*ae_k)+(jc_h[1]*ae_h)+(jc_g[1]*ae_g))            
      )
      
      
      #LN
      O_S_conc[1]=ifelse((jo_k[1]+jo_h[1]+jo_g[1])==0,O_S1[1],  #concentration before inspection
                         O_S[1]/((jo_k[1]*ae_k)+(jo_h[1]*ae_h)+(jo_g[1]*ae_g))
      )
      
      O_conc[1]=ifelse((jo_k[1]+jo_h[1]+jo_g[1])==0,O_S1[1],  #concentration after inspection
                       O[1]/((jo_k[1]*ae_k)+(jo_h[1]*ae_h)+(jo_g[1]*ae_g))            
      )
      
      ##Prevalence
      #Carcass
      prev_CS[1]=ifelse(C_S_conc[1]>0,1,0) #Before inspection
      
      prev_C[1]=ifelse(C_conc[1]>0,1,0) #after inspection
      
      #LN
      prev_OS[1]=ifelse(O_S_conc[1]>0,1,0) #Before inspection
      
      prev_O[1]=ifelse(O_conc[1]>0,1,0) #after inspection
      
      
      ##logcfu/cm^2
      #Carcass
      logCS[1]=ifelse(C_S_conc[1]==0,NA, log(C_S_conc[1])) #Before inspection
      logC[1]=ifelse(C_conc[1]==0,NA, log(C_conc[1])) #after inspection
      
      #LN
      logOS[1]=ifelse(O_S_conc[1]==0,NA,log(O_S_conc[1])) #Before inspection
      logO[1]=ifelse(O_conc[1]==0,NA,log(O_conc[1])) #after inspection
      
      
      ################################
      ##Begin of the carcass looping##
      ################################
      
      for (i in 2:ncarc){  
        
        ## Environment
        
        j_k[i]=round(sample(co,size=1,prob=p_k,replace=T),0) #total number of cuts
        j_h[i]=round(sample(co,size=1,prob=p_h,replace=T),0) #total number of handling
        j_g[i]=round(sample(co,size=1,prob=p_g,replace=T),0) #total number of hook
        
        
        jo_k[i]=rbinom(1,j_k[i],0.5) #number of cuts organs
        jo_h[i]=rbinom(1,j_h[i],0.5) #number of handling organs
        jo_g[i]=rbinom(1,j_g[i],0.5) #number of hooks organs
        
        
        jc_k[i]=j_k[i]-jo_k[i] #number of cuts carcass
        jc_h[i]=j_h[i]-jo_h[i] #number of handling carcass
        jc_g[i]=j_g[i]-jo_g[i] #number of hook carcass
        
        
        
        ##Carcass
        ##S-1
        C_S1[i]=rpois(1,10^rnorm(1,mu[[k]][w],sig[[k]][w]))*rbinom(1,1,prev) #carcass concentration surface
        
        C_kS[i]=(rpois(1,C_S1[i]*ae_k*jc_k[i])) #cfu in knife area
        C_hS[i]=(rpois(1,C_S1[i]*ae_h*jc_h[i])) #cfu in hand area
        C_gS[i]=(rpois(1,C_S1[i]*ae_g*jc_g[i])) #cfu in hook area
        
        ##S
        #Detection
        C_k1[i]=rbinom(1,C_kS[i],1-dc) #cfu after detection knife area
        C_h1[i]=rbinom(1,C_hS[i],1-dc) #cfu after detection hand area
        C_g1[i]=rbinom(1,C_gS[i],1-dc) #cfu after detection hook area
        
        #Loose to environment
        C_k[i]=rbinom(1,C_k1[i],(1-bc_k)^jc_k[i]) #cfu 1-losse to knife
        C_h[i]=rbinom(1,C_h1[i],(1-bc_h)^jc_h[i]) #cfu 1-loose to hand
        C_g[i]=rbinom(1,C_g1[i],(1-bc_g)^jc_g[i]) #cfu 1-loose to hook
        
        
        ## Organ
        ##S-1
        O_S1[i]=rpois(1,rtriangle(n=1,a=min,b=max[[k]][w],c=mod))*rbinom(1,1,alpha_o[[k]][w]/(alpha_o[[k]][w]+beta_o)) # LN concentration
        
        O_kS[i]=rpois(1,O_S1[i]*ae_k*jo_k[i]) #cfu in knife area
        O_hS[i]=rpois(1,O_S1[i]*ae_h*jo_h[i]) #cfu in hand area
        O_gS[i]=rpois(1,O_S1[i]*ae_g*jo_g[i]) #cfu in hook area
        
        ##S
        #Detection
        O_k1[i]=rbinom(1,O_kS[i],1-do)  #cfu after detection knife area
        O_h1[i]=rbinom(1,O_hS[i],1-do)  #cfu after detection hand area
        O_g1[i]=rbinom(1,O_gS[i],1-do)  #cfu after detection hook area
        
        #Loose environment
        O_k[i]=rbinom(1,O_k1[i],(1-bo_k)^jo_k) #cfu 1-losse to knife
        O_h[i]=rbinom(1,O_h1[i],(1-bo_h[[k]][w])^jo_h[i]) #cfu 1-losse to hand
        O_g[i]=rbinom(1,O_g1[i],(1-bo_g)^jo_g[i]) #cfu 1-losse to hook
        
        
        ##Environment
        
        #Knife
        if (check.integer(i/24)=='TRUE'){
          K[i]=(K[i-1]-(O_k2[i-1]+C_k2[i-1])+(C_k1[i]-C_k[i])+(O_k1[i]-O_k[i]))*rbinom(1,1,1-ck)
        }else{K[i]=K[i-1]-(O_k2[i-1]+C_k2[i-1])+(C_k1[i]-C_k[i])+(O_k1[i]-O_k[i])
        }
        
        #Hands
        H[i]=H[i-1]-(O_h2[i-1]+C_h2[i-1])+(C_h1[i]-C_h[i])+(O_h1[i]-O_h[i])
        
        #Hook
        G[i]=G[i-1]-(O_g2[i-1]+C_g2[i-1])+(C_g1[i]-C_g[i])+(O_g1[i]-O_g[i])
        
        
        ##Carcass  
        #Gain from environment
        C_k2[i]= K[i]-rbinom(1,K[i],(1-ak_c)^jc_k[i]) #gain form knife
        C_h2[i]= H[i]-rbinom(1,H[i],(1-ah_c[[k]][w])^jc_h[i]) #gain form hands 
        C_g2[i]= G[i]-rbinom(1,G[i],(1-ag_c)^jc_g[i]) #gain form hook
        
        C_S[i]=C_kS[i]+C_hS[i]+C_gS[i] #Surface before inspection
        C[i]=C_k[i]+C_k2[i]+C_h[i]+C_h2[i]+C_g[i]+C_g2[i] #Surface after insepction
        
        
        ##Organ
        #Gain from environment
        O_k2[i]= (K[i]-C_k2[i])-rbinom(1,(K[i]-C_k2[i]),(1-ak_o)^jo_k[i]) #gain form knife
        O_h2[i]= (H[i]-C_h2[i])-rbinom(1,(H[i]-C_h2[i]),(1-ah_o[[k]][w])^jo_h[i]) #gain form hands 
        O_g2[i]= (G[i]-C_g2[i])-rbinom(1,(G[i]-C_g2[i]),(1-ag_o)^jo_g[i]) #gain form hook
        
        O_S[i]=O_kS[i]+O_hS[i]+O_gS[i] #LN before inspection
        O[i]=O_k[i]+O_k2[i]+O_h[i]+O_h2[i]+O_g[i]+O_g2[i] #LN after inspection
        
        
        
        ##Outputs
        
        
        ##cfu/cm^2
        #Crcass
        C_S_conc[i]=ifelse((jc_k[i]+jc_h[i]+jc_g[i])==0,C_S1[i],  #concentration before inspection
                           C_S[i]/((jc_k[i]*ae_k)+(jc_h[i]*ae_h)+(jc_g[i]*ae_g))
        )
        
        C_conc[i]=ifelse((jo_k[i]+jo_h[i]+jo_g[i])==0,C_S1[i],  #concentration after inspection
                         C[i]/((jc_k[i]*ae_k)+(jc_h[i]*ae_h)+(jc_g[i]*ae_g))            
        )
        
        
        
        #LN
        O_S_conc[i]=ifelse((jo_k[i]+jo_h[i]+jo_g[i])==0,O_S1[i],  #concentration before inspection
                           O_S[i]/((jo_k[i]*ae_k)+(jo_h[i]*ae_h)+(jo_g[i]*ae_g))
        )
        
        O_conc[i]=ifelse((jo_k[i]+jo_h[i]+jo_g[i])==0,O_S1[i],  #concentration after inspection
                         O[i]/((jo_k[i]*ae_k)+(jo_h[i]*ae_h)+(jo_g[i]*ae_g))            
        )
        
        
        
        ##Prevalence
        #Carcass
        prev_CS[i]=ifelse(C_S_conc[i]>0,1,0) #Before inspection
        
        prev_C[i]=ifelse(C_conc[i]>0,1,0) #after inspection
        
        #LN
        prev_OS[i]=ifelse(O_S_conc[i]>0,1,0) #Before inspection
        
        prev_O[i]=ifelse(O_conc[i]>0,1,0) #after inspection
        
        
        ##logcfu/cm^2
        #Carcass
        logCS[i]=ifelse(C_S_conc[i]==0,NA, log(C_S_conc[i])) #Before inspection
        logC[i]=ifelse(C_conc[i]==0,NA, log(C_conc[i])) #after inspection
        
        #LN
        logOS[i]=ifelse(O_S_conc[i]==0,NA,log(O_S_conc[i])) #Before inspection
        logO[i]=ifelse(O_conc[i]==0,NA,log(O_conc[i])) #after inspection
        
        
      } # End of the carcass loop i
      
      me_CS[j]=ifelse(is.na(mean(logCS,na.rm=T)),0,mean(logCS,na.rm=T))
      me_C[j]=ifelse(is.na(mean(logC,na.rm=T)),0,mean(logC,na.rm=T))
      me_OS[j]=ifelse(is.na(mean(logOS,na.rm=T)),0,mean(logOS,na.rm=T))
      me_O[j]=ifelse(is.na(mean(logO,na.rm=T)),0,mean(logO,na.rm=T))
      
      prevalence_CS[j]=mean(prev_CS,na.rm=T)
      prevalence_C[j]=mean(prev_C,na.rm=T)
      prevalence_OS[j]=mean(prev_OS,na.rm=T)
      prevalence_O[j]=mean(prev_O,na.rm=T)
      
      
    } # End of the stochastic node j
    
    mu_mu_c[[k]][w]<-mean(me_C,na.rm=T)
    mu_p_c[[k]][w]<-mean(prevalence_C,na.rm=T)
      
        } #End of the loop in w (internal list loop)
  } #End of scenario looping K (external list loop)
  
  names(mu_mu_c)<-names[1:k]
  names(mu_p_c)<-names[1:k]
  

  dir.create(here("Cap_1","Output","Uni_scenarios"),showWarnings = F)
  
  
  data.local1<- file.path(paste(here(),"/Cap_1","/Output","/Uni_scenarios",sep=""),"mu_mu")
  data.local2<- file.path(paste(here(),"/Cap_1","/Output","/Uni_scenarios",sep=""),"mu_p")

  
  capture.output(mu_mu_c, file = data.local1, append=FALSE)
  capture.output(mu_p_c, file = data.local2, append=FALSE)
  

  
  dir.create(here("Cap_1","Figures","Uni_scenarios"),showWarnings = F)
  
  
  data1<-cbind.data.frame(y=unlist(mu_mu_c[1:4]),x=rep(x,4),
                          parameter=c(
                          rep("Mean[Ci,(S-1)]",length(x)),
                          rep("Standard deviation[Ci,(S-1)]",length(x)),
                          rep("Prev Oi,(S-1)",length(x)),  
                          rep("Oi,C(S-1) Max",length(x))
                          )
                          
  )
  
  #png(file=here("Cap_1","Figures", "Uni_scenarios",paste("Fig3",".png",sep="")),
  #    width = 465, height = 225, units='mm', res = 300)
  
  
  p1<-ggplot(data1,aes(x=x, y=y,color=parameter,shape=parameter))+
    geom_line(size=1.2)+
    geom_point(size=3)+
    scale_colour_grey()+
    scale_shape_manual(values=c(15, 17,18,16))+
    theme_classic()+
    xlab(expression(paste(italic("x"), "value")))+
    ylab(expression(paste("Mean of means ",mu," on contaminaminated carcasses surface "," logCFU","/",cm^2,")")))+
    scale_x_continuous(breaks=seq(-1,1,0.2))+
    theme(panel.grid = element_blank(),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text.x=element_text(size=rel(2)),
          axis.text.y=element_text(size=rel(2)),
          legend.title = element_blank(),
          legend.position="bottom",
          legend.text = element_text(size=rel(1.2)))

  ggsave(filename=here("Cap_1","Figures", "Uni_scenarios","Fig3.png"),p1,
              dpi = 300, height = 10, width = 15, unit = 'in',device="png",limitsize = F)  
  
                          
  data2<-cbind.data.frame(y=unlist(mu_mu_c[5:7]),x=rep(x,3),
                          parameter=c(rep("b O,H",length(x)),
                          rep("a H,O",length(x)),
                          rep("a H,C",length(x))  )
  )
  
  
  p2<-ggplot(data2,aes(x=x, y=y,color=parameter,shape=parameter))+
    geom_line(size=1.2)+
    geom_point(size=3)+
    scale_colour_grey()+
    scale_shape_manual(values=c(15, 17,16))+
    theme_classic()+
    xlab(expression(paste(italic("x"), "value")))+
    ylab(expression(paste("Mean of means ",mu," on contaminaminated carcasses surface "," logCFU","/",cm^2,")")))+
    scale_x_continuous(breaks=seq(-1,1,0.2))+
    theme(panel.grid = element_blank(),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text.x=element_text(size=rel(2)),
          axis.text.y=element_text(size=rel(2)),
          legend.title = element_blank(),
          legend.position="bottom",
          legend.text = element_text(size=rel(1.2)))
  ggsave(filename=here("Cap_1","Figures", "Uni_scenarios","Fig4.png"),p2,
         dpi = 300, height = 10, width = 15, unit = 'in',device="png",limitsize = F)  
  
    
  data3<-cbind.data.frame(y=unlist(mu_p_c[c(3,5:7)]),x=rep(x,4),
                          parameter=c(rep("Prev Oi,(S-1)",length(x)), 
                          rep("b O,H",length(x)),
                          rep("a H,O",length(x)),
                          rep("a H,C",length(x))  )
  )
  
  
  
  p3<-ggplot(data3,aes(x=x, y=y,color=parameter,shape=parameter))+
    geom_line(size=1.2)+
    geom_point(size=3)+
    scale_colour_grey()+
    scale_shape_manual(values=c(15, 17, 18,16))+
    theme_classic()+
    xlab(expression(paste(italic("x"), "value")))+
    ylab("Prevalence of contaminated carcass surface")+
    scale_x_continuous(breaks=seq(-1,1,0.2))+
    theme(panel.grid = element_blank(),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text.x=element_text(size=rel(2)),
          axis.text.y=element_text(size=rel(2)),
          legend.title = element_blank(),
          legend.position="bottom",
          legend.text = element_text(size=rel(1.2)))
  
  ggsave(filename=here("Cap_1","Figures", "Uni_scenarios","Fig5.png"),p3,
         dpi = 300, height = 10, width = 15, unit = 'in',device="png",limitsize = F)  
  
  
  } #End of function

