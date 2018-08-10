NEURON
{
  SUFFIX naIn 
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

  ah = -0.11096151393643595     (/mV) 
  bh = 5.431610352650128     (1) 
  vhh = -102.05330509235864     (mV) 
  Ah = 2.1354326537972503     (/ms) 
  b1h = -0.04414413439854798     (/mV) 
  c1h = 0.0007893023139076368     (/mV2) 
  d1h = -4.4002387790870165e-06     (/mV3) 
  b2h = 0.044148711073748646     (/mV) 
  c2h = -0.007768330848679123     (/mV2) 
  d2h = 6.674307742794694e-05     (/mV3) 

  am = 0.21739294422962296     (/mV) 
  bm = -11.43491119902939     (1) 
  vhm = -99.72043708780117     (mV) 
  Am = 2.0197673427077176     (/ms) 
  b1m = 0.0023509198073411717     (/mV) 
  c1m = 1.650209215350798e-05     (/mV2) 
  d1m = 3.589398825184556e-08     (/mV3) 
  b2m = 0.0023210371050760114     (/mV) 
  c2m = 1.1933380911627887e-05     (/mV2) 
  d2m = -4.72660533814294e-08     (/mV3) 
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
  g = gbar*h*m
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