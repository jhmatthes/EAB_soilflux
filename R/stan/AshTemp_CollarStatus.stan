/* Test a Stan mixed effects model with the ash plot data
*/

data {
  int<lower=1> N;                   //number of data points
  real temp[N];                       //soil moisture
  int<lower=1> C;                   //number of collars
  int<lower=1> S;                   //number of statuses
  int<lower=1> M;                   //number of months
  int<lower=1, upper=C> collar[N];  //collar id
  int<lower=1, upper=M> month[N];   //month id
  int<lower=1, upper=S> status[N];  //status id
}

parameters {
  vector[C] collar_int;      //collar intercepts
  vector[M] month_int;       //month intercepts
  vector[S] status_int;      //status intercepts
  real<lower=0> sigma_e;     //error sd
  real<lower=0> sigma_c;     //collar sd
  real<lower=0> sigma_s;     //status sd
  real<lower=0> sigma_m;     //status sd
}

model {
  real mu;
  //priors
  collar_int ~ normal(0,sigma_c);   //collar random effects
  status_int ~ normal(0,sigma_s);    //status random effects
  month_int ~ normal(0,sigma_m);    //status random effects

  // likelihood
  for (i in 1:N){
    mu = collar_int[collar[i]] + status_int[status[i]] + month_int[month[i]];
    temp[i] ~ normal(mu,sigma_e);
  }
}
