data {
  int<lower=0> N; // number of respondents
  int<lower=0> M; // number of items
  int<lower=0> C; // number of ordinal responses
  int<lower=0> Nstar; // number of test theta
  vector[Nstar] xstar; // test thetas
  int<lower=0, upper=C> y[N, M]; // N by M response matrix with NA 
}

parameters {
  vector[N] theta; // respondent latent ability
  vector[M] discrimination; // item discrimination parameters
  ordered[C-1] difficulty[M]; // item difficulty parameters
  
}

transformed parameters{
   real<lower=0,upper=1> pstar[Nstar,M,(C-1)]; // latent pstar
   real<lower=0,upper=1> p[Nstar,M,C]; // latent p
   real icc[Nstar, M]; // icc
   for(i in 1:Nstar){
    for(j in 1:M){
      for(c in 1:(C-1)){
        pstar[i,j,c] = inv_logit(difficulty[j,c] - discrimination[j]*xstar[i]);
      }
    }
  }

 for(i in 1:Nstar){
    for(j in 1:M){
      p[i,j,1] = pstar[i,j,1];
      for(c in 2:(C-1)){
        p[i,j,c] = pstar[i,j,c] - pstar[i,j,(c-1)];
      }
      p[i,j,C] = 1 - pstar[i,j,(C-1)];
    }
  }
  
  for(i in 1:Nstar){
    for(j in 1:M){
      icc[i,j] = 0;
      for(c in 1:C){
        icc[i,j] += p[i,j,c]*c;
      }
    }
  }
}

model {
  theta ~ normal(0, 1);
  discrimination ~ normal(0, 0.5);
  for(j in 1:M){
    difficulty[j] ~ normal(0, 0.5);
  }
  for(i in 1:N){
    for(j in 1:M){
      if(y[i,j]>0){
        y[i,j] ~ ordered_logistic(theta[i]*discrimination[j], difficulty[j]);
      }
    }
  }
}

