data {  
  int N;
  int streetId[N];
  int otherStreetIndex[N];
  int count[N]; 
}
parameters {
  real intensity[N];
  real<lower=0> stdDeviation;
}
model {
  stdDeviation ~ lognormal(0,1);
  for (i in 1:N) {
    if (streetId[i] < otherStreetIndex[i]) // avoid overgenerating the data
      count[i] ~ poisson(exp(intensity[1+streetId[i]]) + exp(intensity[1+otherStreetIndex[i]]));
    if (otherStreetIndex[i] == -1)
      count[i] ~ poisson(exp(intensity[1+streetId[i]]));
    if (i < N)
      intensity[i+1] ~ normal(intensity[i], stdDeviation);
  }
}