library(EpiEstim);

geometric<-function(x) exp(sum(log(x))/length(x));

tasasr <- function(M, m1, m2, Ncountry, mu, gammap, name_country, tipperpoint, mean_infect, desv_infect)
{ 
  alpha <- NULL;
  G <- NULL;
  R0 <- NULL;
  R02 <- NULL;
  
  Scountry = Ncountry - (M[,3]) - (M[,2]) - (M[,1]);
  NP_country =Ncountry-(M[,3]);
  P_country = M/NP_country;
  Sf_country = Scountry/NP_country;
  m=(m1:m2);
  nm=m2-m1+1;
  
  nch = dim(M);
  
  Ncase <-matrix(0,(nch[1]),1);
  Ncase[1]<-M[1,1];
  for (i in 2:(nch[1])){
    Ncase[i]=M[i,1]-M[i-1,1];
  }
  colnames(Ncase) <- c("New_cases");
  rownames(Ncase) <- rownames(M)[1:(nch[1])];
  #print(Ncase);
  #plot(Ncase);
  dates <- as.data.frame(rownames(M)[1:(nch[1])]);
  rownames(dates) <- rownames(M)[1:(nch[1])];
  
  #Ncase <- as.data.frame(Ncase);
  R0_R0 <- matrix(0,(nch[1]-1),2);
  res <- estimate_R(as.data.frame(Ncase), method = "parametric_si", config = make_config(list(mean_si = mean_infect, std_si = desv_infect)));
  R0_R0[res[[1]]$t_end-1,1] <- res[[1]]$`Mean(R)`;
  R0_R0[res[[1]]$t_end-1,2] <- res[[1]]$`Std(R)`;
  rownames(R0_R0) <- rownames(M)[1:(nch[1]-1)];
  colnames(R0_R0) <- c(paste(name_country,"_Mean_Rt_cori_meth",sep=""),paste(name_country,"_Std_Rt_cori_meth",sep=""));
  
  method1 <- matrix(0,(nch[1]),2);
  for(i in 1:(nch[1])){
    #method1[i,1]=M[i,1]+Scountry[i]-Ncountry+M[i,3];
    method1[i,2]=(NP_country[i]-Scountry[i]-M[i,1])/NP_country[i];
    method1[i,1]=log(Scountry[i]/Scountry[1]);
    #method1[i,2]=log(NP_country[i]/Ncountry);
  }
  colnames(method1) <- c("SI","logS");
  #print(method1);
 
  #method1 <- matrix(0,(nch[1]-1),3);
  #for(i in 1:(nch[1]-1)){
  #  method1[i,1]=M[i+1,1];
  #  method1[i,2]=M[i,1];
  #  method1[i,3]=(M[i,1])^2;
  #}
  #colnames(method1) <- c("x_n+1","x_n","x2_n");
  #plot(method1[,2],method1[,1]);
  #abline(a=0, b=(method1[nch[1],1]/method1[nch[1],2]));

  alpha_brauer <- matrix(0,(nch[1]-1),1);
  R0_brauer <- matrix(0,(nch[1]-1),1);
  Reff_moment <- matrix(0,(nch[1]-1),1);
  for(i in 2:(nch[1])){
    #modelrm <- lm(x_n+1 ~ x_n+x2_n-1, data=as.data.frame(method1[1:i,]), );
    if (i<=tipperpoint){
      data1=as.data.frame(method1[1:i,]);
    } else {
      data1=as.data.frame(method1[tipperpoint:i,]);
    }
    #print(data1);
    modelrm <- lm(SI ~ logS-1, data=data1);
    #print(summary(modelrm));
    alpha_brauer[i-1]=-1*(gammap)*(modelrm$coefficients[1]);
    if(is.na(alpha_brauer[i-1])){
      alpha_brauer[i-1] <- 0;
    }
    if(is.infinite(alpha_brauer[i-1])){
      alpha_brauer[i-1] <- 0;
    }
    R0_brauer[i-1]=(alpha_brauer[i-1])/(gammap + mu);
    Reff_moment[i-1]=(Ncase[i]/(gammap*M[i,1]));
  }
  alpha_brauer = as.data.frame(alpha_brauer);
  colnames(alpha_brauer) <- c(paste(name_country,"_beta_weiss_meth",sep=""));
  rownames(alpha_brauer) <- rownames(M)[1:(nch[1]-1)];
  R0_brauer = as.data.frame(R0_brauer);
  colnames(R0_brauer) <- c(paste(name_country,"_R0_beta_weiss_meth",sep=""));
  rownames(R0_brauer) <- rownames(M)[1:(nch[1]-1)];
  Reff_moment = as.data.frame(Reff_moment);
  colnames(Reff_moment) <- c(paste(name_country,"_R_effec_crespo_meth",sep=""));
  rownames(Reff_moment) <- rownames(M)[1:(nch[1]-1)];
  
  finalc <- cbind(dates[1:(nch[1]-1),],M[1:(nch[1]-1),],Ncase[1:(nch[1]-1),],alpha_brauer,R0_brauer,Reff_moment,R0_R0);
  colnames(finalc)[1] <- c("dates");
  colnames(finalc)[5] <- c("New_cases_confirmed");
  
  alpha <- matrix(0,(nch[1]-1),1);
  R0 <- matrix(0,(nch[1]-1),1);
  for(i in 1:(nch[1]-1)){
    alpha[i]=(NP_country[1]*(Scountry[i]-Scountry[i+1]))/(Scountry[i]*M[i,1]);
    R0[i]=(alpha[i])/(gammap + mu);
  }
  cont_edo=alpha[nch[1]-1]*Sf_country[nch[1]-1]*P_country[nch[1]-1,1]*NP_country[nch[1]-1];
  alpha = as.data.frame(alpha);
  colnames(alpha) <- c(paste(name_country,"_beta_edo",sep=""));
  rownames(alpha) <- rownames(M)[1:(nch[1]-1)];
  R0 = as.data.frame(R0);
  colnames(R0) <- c(paste(name_country,"_R0_beta_edo",sep=""));
  rownames(R0) <- rownames(M)[1:(nch[1]-1)];
  
  colnames(M) <- paste(name_country,colnames(M),sep="_");
  
  finalc <- cbind(finalc,alpha,R0);
  
  #inew <- NULL;
  #inew[1] <- M[1,1];
  #for(i in 2:(nch[1])){
  #  inew[i]=M[i,1]-M[i-1,1];
  #}

  l<- NULL;
  l[1] <-0;
  for(i in 2:(nch[1]-1)){
    l[i]=(((M[i,1])/(M[i-1,1]))-1);
    #if(is.na(l[i])){
    #  l[i] <- 0;
    #}
    #if(is.infinite(l[i])){
    #  l[i] <- 0;
    #}
  }
  l1<- NULL;
  l1[1] <-0;
  for(i in 2:(nch[1]-1)){
    l1[i]=(((M[i,1])/(M[i-1,1])));
  }

  l=as.data.frame(l);
  l1=as.data.frame(l1);
  nl=dim(l);
  
  G <- matrix(0,(nch[1]-1),nm);
  Gg <- matrix(0,(nch[1]-1),nm);
  R02 <- matrix(0,(nch[1]-1),nm);
  Rg02 <- matrix(0,(nch[1]-1),nm);
  for (j in m){
    for(i in 2:(nl[1]-j+1)){
      G[i+j-1,j-m1+1]=mean(l[i:(i+(j-1)),]);
      aux=geometric(l1[i:(i+(j-1)),]);
      if(aux!=0){
        Gg[i+j-1,j-m1+1]=aux-1;
      }
    }
  }
  colnames(G) <- c(paste(name_country,"Gd_m",as.character(m),sep="_"));
  rownames(G) <- rownames(M)[1:(nch[1]-1)];
  G <- as.data.frame(G);
  colnames(Gg) <- c(paste(name_country,"Gdg_m",as.character(m),sep="_"));
  rownames(Gg) <- rownames(M)[1:(nch[1]-1)];
  Gg <- as.data.frame(Gg);
  
  for(i in 1:(nm)){
    R02[,i]=(G[,i])/(gammap + mu);
    Rg02[,i]=(Gg[,i])/(gammap + mu);
  }
  colnames(R02) <-c(paste(name_country,"R0_g_m",as.character(m),sep="_"));
  rownames(R02) <- rownames(M)[1:(nch[1]-1)];
  colnames(Rg02) <- c(paste(name_country,"R0_Gg_m",as.character(m),sep="_"));
  rownames(Rg02) <- rownames(M)[1:(nch[1]-1)];
  R02 = as.data.frame(R02);
  Rg02 = as.data.frame(Rg02);
  
  a<-as.matrix(M[nch[1]-1,1]*(1+G[nch[1]-1,]))-M[nch[1],1];
  b<-as.matrix(M[nch[1]-1,1]*(1+Gg[nch[1]-1,]))-M[nch[1],1];
  
  a1 <- which(max(a)==a);
  b1 <- which(max(b)==b);
  
  Reff <- R0_brauer;
  for(i in 1:(nch[1]-1)){
    Reff[i,]=(Reff[i,])*(Scountry[i])/NP_country[i];
  }
  Reff=as.data.frame(Reff);
  colnames(Reff) <- c("R_effective");
  
  finalc <- cbind(finalc,Reff);
  
  finalc <- as.data.frame(cbind(finalc,G[,a1], R02[,a1], Gg[,b1], Rg02[,b1]));
  nf <- dim(finalc);
  colnames(finalc)[(nf[2]-3):nf[2]] <- cbind(colnames(G)[a1],colnames(R02)[a1],colnames(Gg)[b1],colnames(Rg02)[b1])
  
  return(finalc);
};
