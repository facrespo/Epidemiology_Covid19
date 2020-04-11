library(deSolve);

epsilon=1;


mu= 0.016; #Born rate
#mu= 0;
N = 19000000;
alpha=0.19; #Transmission rate
#alpha=0.3;
beta = (1/3.5);
tmorth=0.0608; # dead natural rate of Chile 
tmcovid=0.0000035; #Dead Rate of COVID
tt=240; #Total time to simulate
tasahospital=0.0001; #Hospitalized

camas_uci=3500;

p = c(mu, N, alpha,  beta, tmorth, tmcovid);

S=N-epsilon;
I=epsilon;
R=0;
D=0;

start_SIR = c(S, I, R, D);

times = seq(0, 240, by = 1);

R0=(alpha/(beta+mu));

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

resol = ode(y=start_SIR, times=times, func=SIR, parms=p);

t=resol[,"time"];
plot(t,resol[,2],type="l",xlab="time",ylab="People",main="Evolution of Infected with COVID-19");
lines(t,resol[,3],col="red");
lines(t,resol[,4],col="blue");
lines(t,resol[,5],col="orange");
legend(x = "topright", legend = c("Susceptible", "Infected", "Recuperated", "Dead"), fill = c("black", "red", "blue","orange"), title = "Populations")

t=resol[,"time"];
plot(t,resol[,3],type="l",xlab="time",ylab="People",main="Infected with COVID-19");

theta=tmorth+tmcovid;

MaxI=max(resol[,3]);
Tmax=which.max(resol[,3]);
MaxR=max(resol[,4]);
Rmax=which.max(resol[,4]);

Porcentaje_Infectado=(resol[94,3]+resol[94,4])/N;

Hospitalizados=tasahospital*(resol[,3]);

camas_uci_disponible=Hospitalizados-camas_uci;
camas_uci_disponible_i =which(1*(camas_uci_disponible>=0)==1);

victimas_adicionales=max(camas_uci_disponible[camas_uci_disponible_i]);

Dia_d=camas_uci_disponible_i[1];

Muertos_estimados=max(resol[Rmax,5]);

total_resultado=cbind(resol,Hospitalizados);

colnames(total_resultado)[2:5] <- c("Susceptible", "Infected", "Recuperated", "Dead");

which.min(total_resultado[47:100,4])

#plot(total_resultado[1:90,4]);


resol[10,]
resol[29,]
Hospitalizados[29]

Hospitalizados[37]
resol[37,]

resol[Dia_d,]

