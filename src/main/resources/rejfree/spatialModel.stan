data {  
  int N;
  int streetId[N];
  int otherStreetIndex[N];
  real distanceToNextInThisStreet[N]; 
  int count[N]; 
}
parameters {
  real intensity[N];
  real<lower=0> slope; 
  real<lower=0> intercept;
}
model {
  slope ~ lognormal(0,1);
  for (i in 1:N) {
    if (streetId[i] < otherStreetIndex[i]) // avoid overgenerating the data
      count[i] ~ poisson(exp(intensity[1+streetId[i]]) + exp(intensity[1+otherStreetIndex[i]]));
    if (distanceToNextInThisStreet[i] != -1.0)
      intensity[i+1] ~ normal(intensity[i], slope * distanceToNextInThisStreet[i] + 0.0001);
  }
}