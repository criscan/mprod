mprod_fits=function(name,graf_opt,tab_opt){


  #Carga funciones y librerias--------------------
  library(MASS)
  library(readxl)
  library(psych)
  library(openxlsx)
  
  
  #ferror_bdin-----------------------
  ferror_bdin=function(data,parini){
    
    K=exp(parini[1])
    r=exp(parini[2])
    sigma=exp(parini[3])
    rho=exp(parini[4])
    p0=exp(parini[5])
    
    Y=datos_in$datos[,2]
    CPUEdat=datos_in$datos[,3]
    
    n=length(Y)
    B=rep(0,1,n)
    CPUEpred=rep(0,1,n)
    
    B[1]=K*p0
    
    for (t in 2:n)
    { 
      B[t]=max(B[t-1] + r/rho*B[t-1]*(1-(B[t-1]/K)^rho) - Y[t-1],0.1)
    }
    
    id=which(CPUEdat>0)
    
    q=exp(mean(log(CPUEdat[id]/B[id])))
    CPUEpred=q*B
    datos_in$priors[2,]=datos_in$priors[2,]+1e-10
    
    suma=sum((1/sigma*(log(CPUEpred[id])-log(CPUEdat[id])))^2)
    sumaprioris=0.5*sum(((log(c(K,r,sigma,rho,p0))-log(datos_in$priors[1,]))/datos_in$priors[2,])^2)
    
    fun=-n*log(0.5/(sigma*sqrt(6.283185)))+suma+sumaprioris
    
    out=list(suma,sumaprioris)
    
    return(fun)
    
  }
  
  #fdinam_bdin--------------------------
  fdinam_bdin=function(data,parini){
    
    K=exp(parini[1])
    r=exp(parini[2])
    sigma=exp(parini[3])
    rho=exp(parini[4])
    p0=exp(parini[5])
    Y=datos_in$datos[,2]
    CPUEdat=datos_in$datos[,3]
    
    
    n=length(Y)
    Biom=rep(0,1,n)
    CPUEpred=rep(0,1,n)
    
    Biom[1]=K*p0
    
    for (t in 2:n)
    { 
      Biom[t]=max(c(Biom[t-1] + r/rho*Biom[t-1]*(1-(Biom[t-1]/K)^rho) - Y[t-1],0.1))
    }
    
    id=which(CPUEdat>0)
    q=exp(mean(log(CPUEdat[id]/Biom[id])))
    CPUEpred=q*Biom
    G=r/rho*Biom*(1-(Biom/K)^rho)
    Fmort=Y/Biom
    
    salidas=data.frame(Y,CPUEdat,CPUEpred,Biom,G,Fmort) # salidas
    
    
    return(salidas)
    
  }
  
  # Archivo con los datos
  
  datos=as.matrix(read_xlsx(name,sheet=1,col_names = TRUE))
  priors=as.matrix(read_xlsx(name,sheet=2,col_names = TRUE))
  
  # Defino los parámetros y opciones iniciales
  K0=priors[1,1]
  r0=priors[1,2]
  sigma=priors[1,3]
  rho=priors[1,4]
  p0=priors[1,5]
  
  parini=log(c(K0,r0,sigma,rho,p0))
  datos_in=list(datos=datos,priors=priors) # datos 
  pars_fin=optim(par=parini,fn=ferror_bdin, data=datos, method="BFGS",hessian=T)
  
  par=pars_fin$par
  vcorrel=solve(pars_fin$hessian)
  cv_par=sqrt(diag(vcorrel))
  LL=pars_fin$value
  sd_par=cv_par*exp(par)
  
  
  
  #Variables de interes--------------------
  outputs=fdinam_bdin(datos,par)
  K=exp(par[1])
  r=exp(par[2])
  sigma=exp(par[3])
  rho=exp(par[4])
  p0=exp(par[5])
  Bmsy=K*(1/(1+rho))^(1/rho)
  MSY=r/rho*Bmsy*(1-(Bmsy/K)^rho)
  Fmsy=MSY/Bmsy
  
  Biom=outputs$Biom
  Fmort=outputs$Fmort
  CPUEdat=outputs$CPUEdat
  CPUEpred=outputs$CPUEpred
  Y=outputs$Y
  G=outputs$G
  
  
  
  # Incertidumbre
  muestra<-mvrnorm(1000,par,vcorrel)
  par_K=exp(muestra[,1])
  par_r=exp(muestra[,2])
  par_sigma=exp(muestra[,3])
  par_rho=exp(muestra[,4])
  par_p0=exp(muestra[,5])
  
  id=which(outputs$CPUEdat>0)
  resid=log(outputs$CPUEdat[id])-log(outputs$CPUEpred[id])
  resid=resid/sd(resid)
  
  par(mfrow = c(2, 2))
  
  id=which(outputs$CPUEdat>0)
  Yrs=datos_in$datos[,1]
  
  if(graf_opt==T){
    
    #Graficos--------------------------
    plot(Yrs[id],CPUEdat[id],main=paste("Model fit  (K=",round(K,0)," r=",round(r,3),")"),ylim = c(0,max(CPUEdat[id])*1.01),
         ylab="Abundance index",xlab="Year",type="b",pch=20,cex=2)
    lines(Yrs,CPUEpred,col="red",lwd=2)
    suave = smooth.spline(Yrs[id], CPUEdat[id], spar=0.5)
    lines(suave,col="green",lwd=2,lty=1)
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    legend("topright",c("data","model","trend"), col=c("black","red","green"),lty=1,lwd=2,
           bty="n")
    
    plot(Yrs[id],resid,main="Residuals std",ylab="Residual",xlab="Year",pch=20,cex=2)
    lines(Yrs[id],resid,type="h",lwd=2,col="blue")
    abline(h = 0, col = "black")
    
    box()
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    
    plot(CPUEpred[id],resid[id],xlab="Predicted values",ylab="Residuals std",pch=20,cex=2,
         main="Predicted vs residuals")
    abline(h=0)
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    
    ml=lm(log(CPUEpred[id])~log(CPUEdat[id]))
    pvalue=summary(ml)$coefficients[2,4]
    correl=cor(log(CPUEpred[id]),log(CPUEdat[id]))
    plot(CPUEdat[id],CPUEpred[id],main=paste("R2=",round(correl^2,2),"  p-value=",format(pvalue,scientific =T)),
         xlab="Observed",ylab="Predicted",pch=20,cex=2)
    lines(CPUEpred[id],CPUEpred[id],type="l",col="green",lwd=2)
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    
    #Fig2-----------------------------------------------------------------------------
    par(mfrow = c(2, 2))
    plot(Yrs,Biom,ylim = c(0,max(Biom)*1.01),type="l",ylab="Biomass",xlab="Year",
         main="Biomass",pch = 16,cex=1,lwd=2)
    text(min(Yrs)+3,Bmsy*1.2,paste("Bmsy=",round(Bmsy,0)),col="red",cex=1)
    abline(h = Bmsy, col = "red",lty = 2,lwd=2)
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    
    
    plot(Yrs,Y,ylim = c(0,max(Y)*1.01),type="l",ylab="Catch",xlab="Year",
         main="Catch",pch = 16,cex=1,lwd=2)
    lines(Yrs,G,ylim = c(0,max(G)*1.01),type="l", col="green",pch = 16,cex=1,lwd=2)
    legend("topright",c("Catch","G"),col=c("black","green"),lty=1,bty="n",lwd=2)
    text(min(Yrs)+3,MSY*1.2,paste("MSY=",round(MSY,0)),col="red",cex=1)
    abline(h = MSY, col = "red",lty = 2,lwd=2)
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    
    plot(Yrs,Fmort,ylim = c(0,max(Fmort)*1.01),type="l",ylab="Fishing mortality",xlab="Year",
         main="Fishing mortality",pch = 16,cex=1,lwd=2)
    text(min(Yrs)+3,Fmsy*1.2,paste("Fmsy=",round(Fmsy,3)),col="red",cex=1)
    abline(h = Fmsy, col = "red",lty = 2,lwd=2)
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    
    
    x=seq(0,K,K/50)
    y=r/rho*x*(1-(x/K)^rho)
    plot(x,y,type="l",lwd=2,ylab="Catch",xlab="Biomass",
         main="Production", ylim=c(0,1.2*MSY))
    text(Bmsy*1.1,MSY*1.1,paste("MSY=",round(MSY,0)),col="red",cex=1)
    abline(v=Bmsy,col = "red",lty = 2,lwd=2)
    abline(h=MSY,col = "red",lty = 2,lwd=2)
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    
    #Kobe----------------------------------------------------------------
    par(mfrow = c(1, 1))
    nyrs=length(Yrs)
    SPR=Biom/K
    Mort_F=Fmort
    target=Bmsy/K
    
    plot(SPR/target,Mort_F/Fmsy,pch = 16,ylab="F/Fmsy",xlab="B/Bmsy", xlim = c(0,max(SPR/target)), ylim = c(0,max(Mort_F/Fmsy)*1.5), 
         type="l",col="black",lty="dashed",main=paste("B/Bmsy=",round(SPR[nyrs]/target,2),
                                                      " F/Fmsy=",round(Mort_F[nyrs]/Fmsy,2)))
    polygon(c(0,1,1,0),c(0,0,1,1),col="yellow1") #amarillo
    polygon(c(1,1.1*max(SPR/target),1.1*max(SPR/target),1),c(0,0,1,1),col="green") #verde
    polygon(c(1,1.1*max(SPR/target),1.1*max(SPR/target),1),c(1,1,1.5*max(Mort_F/Fmsy),1.5*max(Mort_F/Fmsy)),col="yellow1") #amarillo
    polygon(c(0,1,1,0),c(1,1,1.5*max(Mort_F/Fmsy),1.5*max(Mort_F/Fmsy)),col="tomato1") #rojo
    
    lines(SPR/target,Mort_F/Fmsy,pch = 16, type="b",col="black",lty="dashed",)
    lines(SPR[nyrs]/target,Mort_F[nyrs]/Fmsy,type="p",col="blue",pch = 16,cex=2)
    text(SPR/target*.95,Mort_F/Fmsy,paste(Yrs),cex=0.8)
    
    par(mfrow = c(3, 2))
    hist(par_K,20,xlab="K",main=paste("K=",round(K,0)),col="lightblue")
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    abline(v=K,lwd=2,col="red")
    box()
    hist(par_r,20,xlab="r",main=paste("r=",round(r,2)),col="lightblue")
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    abline(v=r,lwd=2,col="red")
    box()
    hist(par_sigma,20,xlab="sigma",main=paste("sigma=",round(sigma,2)),col="lightblue")
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    abline(v=sigma,lwd=2,col="red")
    box()
    hist(par_rho,20,xlab="rho",main=paste("rho=",round(rho,2)),col="lightblue")
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    abline(v=rho,lwd=2,col="red")
    box()
    hist(par_p0,20,xlab="rho",main=paste("p0=",round(p0,2)),col="lightblue")
    grid(nx = NULL, ny = NULL, lty = 2, col = "gray",lwd = 1)
    abline(v=rho,lwd=2,col="red")
    box()
    
    pairs(par_K~par_r+par_sigma+par_rho+par_p0,col="blue",pch=20,cex=0.1,main="Multiple correlations")
    
  }
  
  #Variables de salida----------------
  parametros=data.frame(Bmsy,MSY,Fmsy,Bmsy_K=Bmsy/exp(par)[1])
  
  variables=data.frame(Yrs=Yrs,CPUEobs=round(CPUEdat,2),CPUEpred=round(CPUEpred,2),Yield=Y,Biomass=round(Biom,0),
                       Fmort=round(Fmort,2),Production=round(G,0),B_Bmsy=round(Biom/Bmsy,2),F_Fmsy=round(Fmort/Fmsy,2))
  LLdat=-length(CPUEdat)*log(0.5/(exp(par[3])*sqrt(6.283185)))+sum((1/exp(par[3])*(log(CPUEpred[id])-log(CPUEdat[id])))^2)
  
  tabpars=round(cbind(c(exp(par),LLdat),c(sd_par,0)),3)
  colnames(tabpars)=c('value','sd')
  rownames(tabpars)=c('K','r','sigma','rho','p0','LogL')
  
  Management=round(parametros,2)
  rownames(Management)="value"
  
  if(tab_opt==T){
    
    wb <- createWorkbook()
    addWorksheet(wb, "Table1_Parameters")
    addWorksheet(wb, "Table2_Variables")
    addWorksheet(wb, "Table3_MSYReferences")
    
    
    writeData(wb, sheet = "Table1_Parameters", x = tabpars, rowNames = T )
    writeData(wb, sheet = "Table2_Variables", x = variables, rowNames = F )
    writeData(wb, sheet = "Table3_MSYReferences", x = Management, rowNames = T )
    
    saveWorkbook(wb,paste("mprodOUT_",name), overwrite = TRUE)
    
  }
  
  
  lista=list(Table1_Parameters=tabpars,Table2_Variables=variables,Table3_MSYReferences=Management)
  return(lista)
  
}



