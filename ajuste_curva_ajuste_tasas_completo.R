setwd("E:/investigacion/covid19");

library(deSolve);

epsilon=1;

confirmados <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv");
recuperados <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv");
muertos <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv");
todos <- read.csv("https://covid.ourworldindata.org/data/who/full_data.csv");

ncli <- 49;
nkri <- 5;
nger <- 10;

chile=t(rbind(as.numeric(confirmados[ncli,]),as.numeric(recuperados[ncli,]),as.numeric(muertos[ncli,])));
nc=dim(chile);
chile=as.data.frame(chile[ncli:nc[1],]);
colnames(chile) <- c("confirmed","recovered","deaths");
rownames(chile) <- colnames(confirmados[49,ncli:nc[1]]);

south_korea=t(rbind(as.numeric(confirmados[144,]),as.numeric(recuperados[144,]),as.numeric(muertos[144,])));
nk=dim(south_korea);
south_korea=as.data.frame(south_korea[nkri:nk[1],]);
colnames(south_korea) <- c("confirmed","recovered","deaths");
rownames(south_korea) <- colnames(confirmados[144,nkri:nk[1]]);

germany=t(rbind(as.numeric(confirmados[121,]),as.numeric(recuperados[121,]),as.numeric(muertos[121,])));
ng=dim(germany);
germany=as.data.frame(germany[nger:ng[1],]);
colnames(germany) <- c("confirmed","recovered","deaths");
rownames(germany) <- colnames(confirmados[121,nger:nk[1]]);


nch = dim(chile);
nkr = dim(south_korea);
ngr = dim(germany)

m1=1;
m2=20;
m=(m1:m2);
nm=m2-m1+1;

N_Chile = 19000000;
mu_chile = 0.01271;
N_Korea = 51446201;
mu_korea = 0.007;
N_Germany = 83019213;
mu_germany = 0.0095;
beta_korea = (1/14);
beta_chile = (1/14);
beta_germany = (1/14);

S_Chile = N_Chile - (chile$deaths) - (chile$confirmed) - (chile$recovered);
S_Korea = N_Korea - (south_korea$deaths) - (south_korea$confirmed) - (south_korea$recovered);
S_Germany = N_Germany - (germany$deaths) - (germany$confirmed) - (germany$recovered);

NP_Chile =N_Chile-(chile$deaths);
NP_Korea =N_Korea-(south_korea$deaths);
NP_Germany =N_Germany-(germany$deaths);

P_Chile = chile/NP_Chile;
P_Korea = south_korea/NP_Korea;
P_Germany = germany/NP_Germany;
Sf_Chile = S_Chile/NP_Chile;
Sf_Korea = S_Korea/NP_Korea;
Sf_Germany = S_Germany/NP_Germany;

alpha_chile <- NULL;
G_chile <- NULL;
R0_chile <- NULL;
R02_chile <- NULL;
alpha_korea <- NULL;
R0_korea <- NULL;
R02_korea <- NULL;
G_korea <-NULL
alpha_germany <- NULL;
R0_germany <- NULL;
R02_germany <- NULL;
G_germany <-NULL

alpha_chile <- matrix(0,(nch[1]-1),1);
R0_chile <- matrix(0,(nch[1]-1),1);
for(i in 1:(nch[1]-1)){
  alpha_chile[i]=(P_Chile[i+1,1])/(P_Chile[i,1]*Sf_Chile[i]);
  R0_chile[i]=(alpha_chile[i])/(beta_chile + mu_chile);
}
cont_chile=alpha_chile[nch[1]-1]*Sf_Chile[nch[1]-1]*P_Chile[nch[1]-1,1]*NP_Chile[nch[1]-1];

alpha_chile = as.data.frame(alpha_chile);
colnames(alpha_chile) <- c("Chile_alpha_edo");
rownames(alpha_chile) <- rownames(chile)[1:(nch[1]-1)];
R0_chile = as.data.frame(R0_chile);
colnames(R0_chile) <- c("Chile_R0_alpha_edo");
rownames(R0_chile) <- rownames(chile)[1:(nch[1]-1)];

colnames(chile) <- paste('Chile',colnames(chile),sep="_");

finalc <- cbind(chile[1:(nch[1]-1),],alpha_chile,R0_chile);

l<- NULL;
for(i in 2:(nch[1]-1)){
  l[i-1]=(((P_Chile[i,1])/(P_Chile[i-1,1]))-1);
}

l=as.data.frame(l);
nl=dim(l);

G_chile <- matrix(0,(nch[1]-1),nm);
R02_chile <- matrix(0,(nch[1]-1),nm);
for (j in m){
   for(i in 1:(nl[1]-j)){
     G_chile[i+j+1,j-m1+1]=mean(l[i:(i+(j-1)),]);
     }
  }
colnames(G_chile) <- paste('Chile_Gd_m',as.character(m),sep="_");
rownames(G_chile) <- rownames(chile)[1:(nch[1]-1)];
G_chile <- as.data.frame(G_chile);

for(i in 1:(nm)){
  R02_chile[,i]=(G_chile[,i])/(beta_chile + mu_chile);
}
colnames(R02_chile) <- paste('Chile_R0_G_m',as.character(m),sep="_");
rownames(R02_chile) <- rownames(chile)[1:(nch[1]-1)];
R02_chile = as.data.frame(R02_chile);

finalc <- as.data.frame(cbind(finalc,G_chile, R02_chile));

t= 1:(nch[1]-1);
plot(t,alpha_chile[,1],type="l",xlab="time",ylab="Indice",main="Transmission rate alpha edo COVID-19 Chile");
lines(t,G_chile[,1],col="red");

plot(R0_chile[,1],type="l",xlab="time",ylab="Indice",main="Evolution R0 with alpha edo COVID-19 Chile");

plot(R02_chile[,1],type="l",xlab="time",ylab="Indice",main="Evolution R0 COVID-19 Gd with m=4 Chile");

alpha_korea <- matrix(0,(nkr[1]-1),1);
R0_korea <- matrix(0,(nkr[1]-1),1);
for(i in 1:(nkr[1]-1)){
  alpha_korea[i]=(P_Korea[i+1,1])/(P_Korea[i,1]*Sf_Korea[i]);
  R0_korea[i]=(alpha_korea[i])/(beta_korea + mu_korea);
}

cont_korea=alpha_korea[nkr[1]-1]*Sf_Korea[nkr[1]-1]*P_Korea[nkr[1]-1,1]*NP_Korea[nkr[1]-1];

alpha_korea = as.data.frame(alpha_korea);
colnames(alpha_korea) <- c("South_Korea_alpha_edo");
rownames(alpha_korea) <- rownames(south_korea)[1:(nkr[1]-1)];
R0_korea = as.data.frame(R0_korea);
colnames(R0_korea) <- c("South_Korea_R0_alpha_edo");
rownames(R0_korea) <- rownames(south_korea)[1:(nkr[1]-1)];

colnames(south_korea) <- paste('South_Korea',colnames(south_korea),sep="_");

finalk <- cbind(south_korea[1:(nkr[1]-1),],alpha_korea,R0_korea);

l<- NULL;
for(i in 2:(nkr[1]-1)){
  l[i-1]=(((P_Korea[i,1])/(P_Korea[i-1,1]))-1);
}

l=as.data.frame(l);
nl=dim(l);

G_korea <- matrix(0,(nkr[1]-1),nm);
R02_korea <- matrix(0,(nkr[1]-1),nm);
for (j in m){
  for(i in 1:(nl[1]-j)){
    G_korea[i+j+1,j-m1+1]=mean(l[i:(i+(j-1)),]);
  }
}
colnames(G_korea) <- paste('South_Korea_Gd_m',as.character(m),sep="_");
rownames(G_korea) <- rownames(south_korea)[1:(nkr[1]-1)];
G_korea <- as.data.frame(G_korea);

for(i in 1:(nm)){
  R02_korea[,i]=(G_korea[,i])/(beta_korea + mu_korea);
}
colnames(R02_korea) <- paste('Korea_R0_G_m',as.character(m),sep="_");
rownames(R02_korea) <- rownames(south_korea)[1:(nkr[1]-1)];
R02_korea = as.data.frame(R02_korea);

finalk <- cbind(finalk,G_korea, R02_korea);

t= 1:(nkr[1]-1);
plot(t,alpha_korea[,1],type="l",xlab="time",ylab="Indice",main="Transmission rate  alpha edo COVID-19 South Korea");
lines(t,G_korea[,1],col="red");

plot(R0_korea[,1],type="l",xlab="time",ylab="Indice",main="Evolution R0  alpha edo COVID-19 South Korea");

plot(R02_korea[,1],type="l",xlab="time",ylab="Indice",main="Evolution R0 COVID-19  Gd with m=4 South Korea");

alpha_germany <- matrix(0,(ngr[1]-1),1);
R0_germany <- matrix(0,(ngr[1]-1),1);
for(i in 1:(ngr[1]-1)){
  alpha_germany[i]=(P_Germany[i+1,1])/(P_Germany[i,1]*Sf_Germany[i]);
  R0_germany[i]=(alpha_germany[i])/(beta_germany + mu_germany);
}

cont_germany=alpha_germany[ngr[1]-1]*Sf_Germany[ngr[1]-1]*P_Germany[ngr[1]-1,1]*NP_Germany[ngr[1]-1];

alpha_germany = as.data.frame(alpha_germany);
colnames(alpha_germany) <- c("Germany_alpha_edo");
rownames(alpha_germany) <- rownames(germany)[1:(ngr[1]-1)];
R0_germany = as.data.frame(R0_germany);
colnames(R0_germany) <- c("Germany_R0_alpha_edo");
rownames(R0_germany) <- rownames(germany)[1:(ngr[1]-1)];

colnames(germany) <- paste('Germany',colnames(germany),sep="_");

finalg <- cbind(germany[1:(ngr[1]-1),],alpha_germany,R0_germany);

l<- NULL;
for(i in 2:(ngr[1]-1)){
  l[i-1]=(((P_Germany[i,1])/(P_Germany[i-1,1]))-1);
}

l=as.data.frame(l);
nl=dim(l);

G_germany <- matrix(0,(ngr[1]-1),nm);
R02_germany <- matrix(0,(ngr[1]-1),nm);
for (j in m){
  for(i in 1:(nl[1]-j)){
    G_germany[i+j+1,j-m1+1]=mean(l[i:(i+(j-1)),]);
  }
}
colnames(G_germany) <- paste('Germany_Gd_m',as.character(m),sep="_");
rownames(G_germany) <- rownames(germany)[1:(ngr[1]-1)];
G_germany <- as.data.frame(G_germany);

for(i in 1:(nm)){
  R02_germany[,i]=(G_germany[,i])/(beta_germany + mu_germany);
}
colnames(R02_germany) <- paste('Germany_R0_G_m',as.character(m),sep="_");
rownames(R02_germany) <- rownames(germany)[1:(ngr[1]-1)];
R02_germany = as.data.frame(R02_germany);

finalg <- cbind(finalg,G_germany, R02_germany);

t= 1:(ngr[1]-1);
plot(t,alpha_germany[,1],type="l",xlab="time",ylab="Indice",main="Transmission rate  alpha edo COVID-19 Germany");
lines(t,G_germany[,1],col="red");

plot(R0_germany[,1],type="l",xlab="time",ylab="Indice",main="Evolution R0 alpha edo COVID-19 Germany");

plot(R02_germany[,1],type="l",xlab="time",ylab="Indice",main="Evolution R0 COVID-19 Gd with m=4 Germany");

nmax=max(nch[1],nkr[1],ngr[1]);

nff = dim(finalc);

tabla_final <- matrix(0,nmax+1,3*nff[2]);
tabla_final[(nmax-nch[1]+1):(nmax),1:nff[2]] <- as.matrix.data.frame(finalc[1:nch[1],1:nff[2]]);
colnames(tabla_final) <- cbind(colnames(finalc),colnames(finalk),colnames(finalg));
tabla_final <- as.data.frame(tabla_final);
rownames(tabla_final)[1:(nmax)] <- c(colnames(confirmados[,nkri:nk[1]]));
rownames(tabla_final)[nmax+1] <- c(paste("Predicted_",rownames(tabla_final)[nmax],sep=""));
tabla_final[nmax,1:3] <- as.matrix(chile[nch[1],1:3]);
tabla_final[nmax+1,4]=cont_chile;
tabla_final[nmax+1,(6):(nm+5)]=as.matrix(chile[nch[1]-1,1]*(1+G_chile[nch[1]-1,]));

tabla_final[(nmax-nkr[1]+1):(nmax),(nff[2]+1):(2*nff[2])] <- as.matrix.data.frame(finalk[1:nkr[1],1:nff[2]]);
tabla_final[nmax,(nff[2]+1):(nff[2]+3)] <- as.matrix(south_korea[nkr[1],1:3]);
tabla_final[nmax+1,(nff[2]+4)]=cont_korea;
tabla_final[nmax+1,(nff[2]+6):(nff[2]+nm+5)]=as.matrix(south_korea[nkr[1]-1,1]*(1+G_korea[nkr[1]-1,]));

tabla_final[(nmax-ngr[1]+1):(nmax),(2*nff[2]+1):(3*nff[2])] <- as.matrix.data.frame(finalg[1:ngr[1],1:nff[2]]);
tabla_final[nmax,(2*nff[2]+1):(2*nff[2]+3)] <- as.matrix(germany[ngr[1],1:3]);
tabla_final[nmax+1,(2*nff[2]+4)]=cont_germany;
tabla_final[nmax+1,(2*nff[2]+6):(2*nff[2]+nm+5)]=as.matrix(germany[ngr[1]-1,1]*(1+G_germany[ngr[1]-1,]));

write.table(tabla_final,file="E:/investigacion/covid19/Final_Table.txt", quote=FALSE, row.names=TRUE, col.names = TRUE, sep=';');

