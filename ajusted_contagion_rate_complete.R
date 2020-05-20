setwd("E:/investigacion/covid19");

library(deSolve);
library(tidyr);
library(dplyr);

source("settingsSIR.R");

confirmados <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv");
recuperados <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv");
muertos <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv");
todos <- read.csv("https://covid.ourworldindata.org/data/who/full_data.csv");

nc <- dim(confirmados);
nr <- dim(recuperados);
nm <- dim(muertos);

confirmados$Country <- as.character(confirmados$Country.Region);
for (i in 1:nc[1]){
  if(nchar(trimws(confirmados$Province.State))[i]!=0){
    confirmados$Country[i] <- (paste(confirmados$Country.Region[i],confirmados$Province.State[i],sep="_"));
  }
}

recuperados$Country <- as.character(recuperados$Country.Region);
for (i in 1:nr[1]){
  if(nchar(trimws(recuperados$Province.State))[i]!=0){
    recuperados$Country[i] <- (paste(recuperados$Country.Region[i],recuperados$Province.State[i],sep="_"));
  }
}

muertos$Country <- as.character(muertos$Country.Region);
for (i in 1:nm[1]){
  if(nchar(trimws(muertos$Province.State))[i]!=0){
    muertos$Country[i] <- (paste(muertos$Country.Region[i],muertos$Province.State[i],sep="_"));
  }
}

country <- c('Chile');

pccl <- which(confirmados$Country==country);
prcl <- which(recuperados$Country==country);
pdcl <- which(muertos$Country==country);

chile=t(rbind(as.numeric(confirmados[pccl,5:nc[2]]),as.numeric(recuperados[prcl,5:nc[2]]),as.numeric(muertos[pdcl,5:nc[2]])));
nc1=dim(chile);
ncli <- which(1*(chile[,1])!=0)[1];
chile=as.data.frame(chile[ncli:nc1[1],]);
colnames(chile) <- c("confirmed","recovered","deaths");
rownames(chile) <- colnames(confirmados[pccl,(ncli+4):nc[2]]);

country2 <- c('Korea, South');
country2.1 <- c('Republic of Korea');

pckr <- which(confirmados$Country==country2);
prkr <- which(recuperados$Country==country2);
pdkr <- which(muertos$Country==country2);

south_korea=t(rbind(as.numeric(confirmados[pckr,5:nc[2]]),as.numeric(recuperados[prkr,5:nc[2]]),as.numeric(muertos[pdkr,5:nc[2]])));
nk=dim(south_korea);
nkri <- which(1*(south_korea[,1])!=0)[1];
south_korea=as.data.frame(south_korea[nkri:nk[1],]);
colnames(south_korea) <- c("confirmed","recovered","deaths");
rownames(south_korea) <- colnames(confirmados[pckr,(nkri+4):nc[2]]);

country3 <- c('Argentina');

pcar <- which(confirmados$Country==country3);
prar <- which(recuperados$Country==country3);
pdar <- which(muertos$Country==country3);

argentina=t(rbind(as.numeric(confirmados[pcar,5:nc[2]]),as.numeric(recuperados[prar,5:nc[2]]),as.numeric(muertos[pdar,5:nc[2]])));
naar=dim(argentina);
nnar <- which(1*(argentina[,1])!=0)[1];
argentina=as.data.frame(argentina[nnar:naar[1],]);
colnames(argentina) <- c("confirmed","recovered","deaths");
rownames(argentina) <- colnames(confirmados[pcar,(nnar+4):nc[2]]);

country4 <- c('Germany');

pcge <- which(confirmados$Country==country4);
prge <- which(recuperados$Country==country4);
pdge <- which(muertos$Country==country4);


germany=t(rbind(as.numeric(confirmados[pcge,5:nc[2]]),as.numeric(recuperados[prge,5:nc[2]]),as.numeric(muertos[pdge,5:nc[2]])));
ng=dim(germany);
nger <- which(1*(germany[,1])!=0)[1];
germany=as.data.frame(germany[nger:ng[1],]);
colnames(germany) <- c("confirmed","recovered","deaths");
rownames(germany) <- colnames(confirmados[pcge,(nger+4):nc[2]]);

country5 <- c('United Kingdom');

pcuk <- which(confirmados$Country==country5);
pruk <- which(recuperados$Country==country5);
pduk <- which(muertos$Country==country5);


uk=t(rbind(as.numeric(confirmados[pcuk,5:nc[2]]),as.numeric(recuperados[pruk,5:nc[2]]),as.numeric(muertos[pduk,5:nc[2]])));
nuk=dim(uk);
nuki <- which(1*(uk[,1])!=0)[1];
uk=as.data.frame(uk[nuki:nc1[1],]);
colnames(uk) <- c("confirmed","recovered","deaths");
rownames(uk) <- colnames(confirmados[pcuk,(nuki+4):nc[2]]);

country6 <- c('Belgium');

pcbe <- which(confirmados$Country==country6);
prbe <- which(recuperados$Country==country6);
pdbe <- which(muertos$Country==country6);


belgium=t(rbind(as.numeric(confirmados[pcbe,5:nc[2]]),as.numeric(recuperados[prbe,5:nc[2]]),as.numeric(muertos[pdbe,5:nc[2]])));
nbe=dim(belgium);
nbei <- which(1*(belgium[,1])!=0)[1];
belgium=as.data.frame(belgium[nbei:nc1[1],]);
colnames(belgium) <- c("confirmed","recovered","deaths");
rownames(belgium) <- colnames(confirmados[pcbe,(nbei+4):nc[2]]);

nch = dim(chile);
nkr = dim(south_korea);
ngr = dim(germany);
nar = dim(argentina);
nuk1 = dim(uk);
nbe1 = dim(belgium);


m1=1;
m2=20;

population <- read.csv("WPP2019_TotalPopulationBySex.csv", sep=",", header=TRUE);
indicators <- read.csv("WPP2019_Period_Indicators_Medium.csv", sep=",", header=TRUE);

N_Chile <- population[(population$Location==country)&(population$VarID==9)&(population$Time==2020),9]*1000;
N_Korea <- population[(population$Location==country2.1)&(population$VarID==9)&(population$Time==2020),9]*1000;
N_Argentina <- population[(population$Location==country3)&(population$VarID==9)&(population$Time==2020),9]*1000;
N_Germany <- population[(population$Location==country4)&(population$VarID==9)&(population$Time==2020),9]*1000;
N_UK <- population[(population$Location==country5)&(population$VarID==9)&(population$Time==2020),9]*1000;
N_Belgium <- population[(population$Location==country6)&(population$VarID==9)&(population$Time==2020),9]*1000;


#(indicators[(indicators$Location==country)&(indicators$VarID==2)&(indicators$MidPeriod==2018),]);

mu_chile <- (indicators[(indicators$Location==country)&(indicators$VarID==2)&(indicators$MidPeriod==2018),9])/1000;
mu_korea <- (indicators[(indicators$Location==country2.1)&(indicators$VarID==2)&(indicators$MidPeriod==2018),9])/1000;
mu_argentina <- (indicators[(indicators$Location==country3)&(indicators$VarID==2)&(indicators$MidPeriod==2018),9])/1000;
mu_germany <- (indicators[(indicators$Location==country4)&(indicators$VarID==2)&(indicators$MidPeriod==2018),9])/1000;
mu_uk <- (indicators[(indicators$Location==country5)&(indicators$VarID==2)&(indicators$MidPeriod==2018),9])/1000;
mu_belgium <- (indicators[(indicators$Location==country6)&(indicators$VarID==2)&(indicators$MidPeriod==2018),9])/1000;

Tau=14;

mean_time <- 7.5 #time mean de presence of sick
desv_time <- 3.4 #desviacion 


gamma_korea = (1/14);
gamma_chile = (1/14);
gamma_germany = (1/14);




chilet <- tasasr(chile, m1, m2, N_Chile, mu_chile, gamma_chile, country, 17, mean_time, desv_time);

write.table(chilet,file="E:/investigacion/covid19/Tasas_Chile_20200519.txt", quote=FALSE, row.names=FALSE, col.names = TRUE, sep=';');

koreat <- tasasr(south_korea, m1, m2, N_Korea, mu_korea, gamma_korea, country2, 1, mean_time, desv_time);

write.table(koreat,file="E:/investigacion/covid19/Tasas_Korea_20200519.txt", quote=FALSE, row.names=FALSE, col.names = TRUE, sep=';');

argentinat <- tasasr(argentina, m1, m2, N_Argentina, mu_argentina, gamma_korea, country3, 1, mean_time, desv_time);

write.table(argentinat,file="E:/investigacion/covid19/Tasas_Argentina_20200519.txt", quote=FALSE, row.names=FALSE, col.names = TRUE, sep=';');

germanyt <- tasasr(germany, m1, m2, N_Germany, mu_germany, gamma_korea, country4, 1, mean_time, desv_time);

write.table(germanyt,file="E:/investigacion/covid19/Tasas_Germany_20200519.txt", quote=FALSE, row.names=FALSE, col.names = TRUE, sep=';');

ukt <- tasasr(uk, m1, m2, N_UK, mu_uk, gamma_korea, country5, 1, mean_time, desv_time);

write.table(germanyt,file="E:/investigacion/covid19/Tasas_UK_20200519.txt", quote=FALSE, row.names=FALSE, col.names = TRUE, sep=';');

belgiumt <- tasasr(belgium, m1, m2, N_Belgium, mu_belgium, gamma_korea, country6, 1, mean_time, desv_time);

write.table(belgiumt,file="E:/investigacion/covid19/Tasas_Belgium_20200519.txt", quote=FALSE, row.names=FALSE, col.names = TRUE, sep=';');


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

write.table(tabla_final,file="E:/investigacion/covid19/Final_Table_20200426.txt", quote=FALSE, row.names=TRUE, col.names = TRUE, sep=';');

