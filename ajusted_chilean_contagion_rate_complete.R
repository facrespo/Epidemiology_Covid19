setwd("E:/investigacion/covid19");

library(xlsx);
library(tidyr);
library(dplyr);

source("settingsSIR.R");

confirmados <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv");
namesc <- colnames(confirmados)[74:121];

new_confirm <- read.xlsx2("contagiados_fern.xlsx", 1, header=TRUE, colClasses=NA);
acum_infected <- read.xlsx2("contagiados_fern.xlsx", 2, header=TRUE, colClasses=NA);
dead <- read.xlsx2("contagiados_fern.xlsx", 3, header=TRUE, colClasses=NA);

ndat <- dim(new_confirm);

colnames(new_confirm)[2:ndat[2]] <- namesc;
colnames(acum_infected)[2:ndat[2]] <- namesc;
colnames(dead)[2:ndat[2]] <- namesc;

hospital <- read.xlsx2("uci_fern.xlsx", 1, header=TRUE, colClasses=NA);
population <- read.xlsx2("uci_fern.xlsx", 2, header=TRUE, colClasses=NA);

colnames(hospital)[2:ndat[2]] <- namesc;

nr <- dim(population);

mu_chile=0.012465;

Tau=14;

m1=1;
m2=20;

mean_time <- 7.5 #time mean de presence of sick
desv_time <- 3.4 #desviacion 

gamma_chile = (1/14);

for (i in 1:nr[1]){
  region <- as.character(population[i,1]);
  N_pol <- population[i,2];
  k1 <- which(new_confirm$regiones==region);
  k2 <- which(acum_infected$regiones==region);
  k3 <- which(dead$regiones==region);
  datat=t(rbind(as.numeric(acum_infected[k1,2:ndat[2]]),0,as.numeric(dead[k3,2:ndat[2]])));
  colnames(datat) <- c("confirmed","recovered","deaths");
  rownames(datat) <- colnames(new_confirm)[2:ndat[2]];
  datatt <- tasasr(datat, m1, m2, N_pol, mu_chile, gamma_chile, region, 1, mean_time, desv_time);
  if(i==1){
    tabla_final <- datatt;
  } else {
    datatt$dates <- NULL;
    tabla_final <- cbind(tabla_final,datatt);
  }
}

write.table(tabla_final,file="E:/investigacion/covid19/Tasas_Region_20200517.txt", quote=FALSE, row.names=FALSE, col.names = TRUE, sep=';');
