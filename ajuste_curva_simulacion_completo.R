setwd("E:/investigacion/covid19");

library(deSolve);

epsilon=1;

confirmados <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv");
recuperados <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv");
muertos <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv");
todos <- read.csv("https://covid.ourworldindata.org/data/who/full_data.csv");


chile=t(rbind(confirmados[49,],recuperados[49,],muertos[49,]));
nc=dim(chile);
chile=chile[46:nc[1],];

korea=t(rbind(confirmados[144,],recuperados[144,],muertos[144,]));
nk=dim(korea);
korea=korea[5:nc[1],];

NSim=5000;

mu= runif(NSim,min=0.014,max=0.017);
N = 19000000;
alpha=runif(NSim,min=0.7,max=0.9);
beta = runif(NSim,min=(1/6),max=(1/3));
tmorth=0.0608;
tmcovid=runif(NSim,min=0.000009,max=0.00004);
tt=240;
tasahospital=runif(NSim,min=0.0005,max=0.003);

camas_uci=3500;

tabla <- matrix(0, NSim, 16);
colnames(tabla) <- c("mu","alpha","beta","rate_death_covid","rate_hospital","R_0","Max_Infected","T_Max_infected","Max_Recupered","T_Max_Recupered","rate_infected","overload_death","Max_hospitalized","T_max_hospitalized","Date_overload","People_death");

SIR = function(t,Z,p){
  S=Z[1]; I=Z[2]; R=Z[3]; D=Z[4]
  N=S+I+R+D;
  mu=p[1]; 
  alpha=p[3]; 
  beta=p[4];
  tmorth=p[5];
  tmcovid=p[6];
  dS=mu*(R)-alpha*S*I/N;
  dI=alpha*S*I/N-(mu+beta)*I;
  dR=beta*I-mu*R;
  dD=tmcovid*I;
  dZ=c(dS,dI,dR,dD);
  return(list(dZ))}

for (i in 1:NSim){

p = c(mu[i], N, alpha[i],  beta[i], tmorth, tmcovid[i]);

tabla[i,1]=mu[i];
tabla[i,2]=alpha[i];
tabla[i,3]=beta[i];
tabla[i,4]=tmcovid[i];
tabla[i,5]=tasahospital[i];

S=N-epsilon;
I=epsilon;
R=0;
D=0;

start_SIR = c(S, I, R, D);

times = seq(0, 240, by = 1);

R0=(alpha[i]/(beta[i]+mu[i]));
tabla[i,6]=R0;

resol = ode(y=start_SIR, times=times, func=SIR, parms=p);

theta=tmorth+tmcovid;

MaxI=max(resol[,3]);
tabla[i,7]=MaxI;
Tmax=which.max(resol[,3]);
tabla[i,8]=Tmax;
MaxR=max(resol[,4]);
tabla[i,9]=MaxR;
Rmax=which.max(resol[,4]);
tabla[i,10]=Rmax;

Porcentaje_Infectado=(resol[Rmax,3]+resol[Rmax,4])/N;
tabla[i,11]=Porcentaje_Infectado;

victimas_adicionales=max(camas_uci_disponible[camas_uci_disponible_i]);
tabla[i,12]=victimas_adicionales;

Hospitalizados <- tasahospital[i]*(resol[,3]);
tabla[i,13]=max(Hospitalizados);
tabla[i,14]=which.max(Hospitalizados);

camas_uci_disponible=Hospitalizados-camas_uci;
camas_uci_disponible_i =which(1*(camas_uci_disponible>=0)==1);

Dia_d=camas_uci_disponible_i[1];
tabla[i,15]=Dia_d;

Muertos_estimados=max(resol[Rmax,5]);
tabla[i,16]=Muertos_estimados;

#total_resultado=cbind(resol,Hospitalizados);

#colnames(total_resultado)[2:5] <- c("Susceptible", "Infected", "Recuperated", "Dead");
}

write.table(tabla,file="Tablafinal.txt",quote=FALSE,row.names=FALSE,sep=';');

