/* Test a Stan mixed effects model with the ash plot data
*/

data {
  int<lower=1> N;                   //number of data points
  real co[N];                       //co2 flux
  real<lower=0,upper=1> sm[N];      //predictor: soil moisture
  real<lower=0,upper=20> temp[N];   //predictor: temperature
  real<lower=0,upper=30> rain[N];   //monthly rainfall
  real<lower=10,upper=100> dbh[N];  //dbh
  real<lower=0,upper=10> root[N];   //root biomass
  real<lower=0.8,upper=1> sand[N];  //fraction sand
  real<lower=1,upper=7> pH[N];      //soil pH
  int<lower=1> C;                   //number of collars
  int<lower=1> M;                   //number of months
  int<lower=1, upper=C> collar[N];  //collar id
  int<lower=1, upper=M> month[N];   //month id
}

parameters {
  vector[8] beta;            //fixed intercept and slopes
  vector[C] collar_int;      //collar intercepts
  vector[M] month_int;       //date intercepts
  real<lower=0> sigma_e;     //error sd
  real<lower=0> sigma_c;     //collar sd
  real<lower=0> sigma_m;     //month sd

}

model {
  real mu;
  //priors
  collar_int ~ normal(0,sigma_c);   //collar random effects
  month_int ~ normal(0,sigma_m);    //month random effects

  // likelihood
  for (i in 1:N){
    mu = beta[1] + collar_int[collar[i]] + month_int[month[i]] + beta[2]*sm[i] + beta[3]*temp[i] + beta[4]*rain[i] + beta[5]*dbh[i] + beta[6]*root[i] + beta[7]*sand[i] + beta[8]*pH[i]; 
    co[i] ~ lognormal(mu,sigma_e);
  }
}
