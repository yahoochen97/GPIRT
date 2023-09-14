// input data
data {
  int<lower=0> n; // num of respondents
  int<lower=0> m; // num of questions
  int<lower=0> horizon; // num of time periods
  int<lower=2> K; // num of ordered categories
  real sigma; // starting dispersion std
  int<lower=0, upper=K> y[n, m, horizon]; // observations
}

// The parameters accepted by the model.
parameters {
  real<lower=-2,upper=2> beta[m, horizon];
  ordered[K-1] alpha[m, horizon];
  real theta[n, horizon];
//  real<lower=0,upper=1> sigma;
}

// The model to be estimated
model {
//  sigma ~ uniform(0,1);
  for(j in 1:m){
      for(h in 1:horizon){
        beta[j,h] ~ normal(0,1);
        for(k in 1:(K-1)){
          alpha[j,h,k] ~ normal(0,2);
        }
      }
  }
  for(i in 1:n){
    theta[i,1] ~ normal(0,1);
    for(h in 2:horizon){
      theta[i,h] ~ normal(theta[i,h-1], sigma);
    }
    for(j in 1:m){
      for(h in 1:horizon){
        if(y[i,j,h]>0){
          y[i,j,h] ~ ordered_logistic(theta[i,h]*beta[j,h], alpha[j,h]);
        }
      }
    }
  }
}

