NPPmodel <- function (rain, temperature, model = "miami") 
{
  if (model == "miami") {
    npprain = 3000 * (1 - exp(-0.000664 * rain))
    npptemp = 3000/(1 + exp(1.315 - 0.119 * temperature))
    npp = pmin(npprain, npptemp)
  }
  else if (model == "schuur") {
    npprain = (0.005215 * (rain)^1.12363)/exp(0.000459532 * 
                                                rain)
    npptemp = 17.6243/(1 + exp(1.3496 - 0.071514 * temperature))
    npp = pmin(npprain, npptemp) * 200
  }
  else if (model == "NCEAS") {
    npprain = (0.551 * (rain)^1.055)/exp(0.000306 * rain)
    npptemp = 2540/(1 + exp(1.584 - 0.0622 * temperature))
    npp = pmin(npprain, npptemp) * 2
  }
  return(npp)
}
