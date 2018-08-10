NEURON
{
  SUFFIX Naf 
  USEION na READ ena WRITE ina 
  RANGE gbar, g, ina
  GLOBAL ena
}

UNITS
{
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)
}

PARAMETER
{
  gbar = 1 (S/cm2)

  ah = -0.1273346355798809     (/mV) 
  bh = 5.152406910303986     (1) 
  vhh = -26.78227119078402     (mV) 
  Ah = 2.7534700687652656     (/ms) 
  b1h = -0.07379058969540574     (/mV) 
  c1h = 0.0016451931990667092     (/mV2) 
  d1h = -8.649261067905068e-06     (/mV3) 
  b2h = -0.0007430238320677331     (/mV) 
  c2h = 0.00039087593418802777     (/mV2) 
  d2h = -2.6476503558852618e-06     (/mV3) 

  am = 0.09667156313267274     (/mV) 
  bm = -1.74740619222468     (1) 
  vhm = -29.47526489831788     (mV) 
  Am = 0.2243602552703792     (/ms) 
  b1m = 0.031965904374381834     (/mV) 
  c1m = 5.8527183893779365e-05     (/mV2) 
  d1m = -3.069187202888886e-06     (/mV3) 
  b2m = 0.013299791036635608     (/mV) 
  c2m = 0.00021221312256497534     (/mV2) 
  d2m = -1.8718728365240889e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m*m*m
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}