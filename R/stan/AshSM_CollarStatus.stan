/* Test a Stan mixed effects model with the ash plot data
*/

data {
  int<lower=1> N;                   //number of data points
  real sm[N];                       //soil moisture
  real<lower=0,upper=30> rain[N];   //monthly rainfall
  real<lower=0,upper=1> sand[N];    //fraction sand
  int<lower=1> C;                   //number of collars
  int<lower=1> S;                   //number of statuses
  int<lower=1, upper=C> collar[N];  //collar id
  int<lower=1, upper=S> status[N];  //status id
}

parameters {
  vector[3] beta;            //fixed intercept and slopes
  vector[C] collar_int;      //collar intercepts
  vector[S] status_int;      //status intercepts
  real<lower=0> sigma_e;     //error sd
  real<lower=0> sigma_c;     //collar sd
  real<lower=0> sigma_s;     //status sd

}

model {
  real mu;
  //priors
  collar_int ~ normal(0,sigma_c);   //collar random effects
  status_int ~ normal(0,sigma_s);    //status random effects

  // likelihood
  for (i in 1:N){
    mu = beta[1] + collar_int[collar[i]] + status_int[status[i]] + beta[2]*rain[i] + beta[3]*sand[i];
    sm[i] ~ normal(mu,sigma_e);
  }
}
