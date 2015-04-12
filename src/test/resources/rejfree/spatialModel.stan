// Terminology: streetAtCorner = ordered pair of streets (the second street might be unk)

data { 
  // Number of latent variables:
  int nStreetAtCorners;
 
  // Geographic information: for each streetAtCorner id, previous streetAtCorner id when moving on the first street in the streetAtCorner pair
  int nPreviousStreetAtCorners;
  int previousStreetAtCornerIndex[nPreviousStreetAtCorners];
  int currentStreetAtCornerIndex[nPreviousStreetAtCorners];
  
  // Accident data
  int nStreetAtCornerAccidents;
  int streetAtCorner1Index[nStreetAtCornerAccidents]; 
  int streetAtCorner2Index[nStreetAtCornerAccidents]; 
  int accidentCount[nStreetAtCornerAccidents];
  
}

parameters {

  real logIntensity[nStreetAtCorners];
  
  // controls spatial smoothing on the streets
  // larger values means that the logIntensity are more free to vary
  real<lower=0> drift;
  
}

model {

  drift ~ lognormal(0,1);
  
  // Geographic prior
  for (geoIndex in 1:nPreviousStreetAtCorners) {
    if (previousStreetAtCornerIndex[geoIndex] == -1) {
      // this is the first street in the series
      logIntensity[1+currentStreetAtCornerIndex[geoIndex]] ~ normal(0,10);
    } else {
      // center to the previous one
      logIntensity[1+currentStreetAtCornerIndex[geoIndex]] ~ normal(logIntensity[1+previousStreetAtCornerIndex[geoIndex]], drift);//0.2);
    }
  }
  
  // Poisson likelihood
  for (accidentIndex in 1:nStreetAtCornerAccidents) {
    if (streetAtCorner2Index[accidentIndex] == -1) {
      // this is a discretization point with zero accidents
      accidentCount[accidentIndex] ~ poisson(exp(logIntensity[1+streetAtCorner1Index[accidentIndex]]));
    } else {
      accidentCount[accidentIndex] ~ poisson(exp(logIntensity[1+streetAtCorner1Index[accidentIndex]]) +
                                             exp(logIntensity[1+streetAtCorner2Index[accidentIndex]]) );
    }
  }

}