NEURON
{
  SUFFIX hh1 
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

  am = 0.10469486949853081     (/mV) 
  bm = -4.143239049233531     (1) 
  vhm = -50.834074632688534     (mV) 
  Am = 0.39034069050329767     (/ms) 
  b1m = 0.027798182726735322     (/mV) 
  c1m = 0.0008199080120512128     (/mV2) 
  d1m = 1.8390737417418868e-05     (/mV3) 
  b2m = 4.36366469269677     (/mV) 
  c2m = -0.0008550992297809958     (/mV2) 
  d2m = -3.67118523108003e-06     (/mV3) 

  ah = -0.14147104998630744     (/mV) 
  bh = 8.793782928508712     (1) 
  vhh = -67.24919895329549     (mV) 
  Ah = 0.605848962392589     (/ms) 
  b1h = 0.0978745815171717     (/mV) 
  c1h = 0.0024844085699179267     (/mV2) 
  d1h = 4.02859779634958e-05     (/mV3) 
  b2h = 0.0902248928298633     (/mV) 
  c2h = -0.0009023835689362094     (/mV2) 
  d2h = 2.7955955298804075e-06     (/mV3) 

  an = 0.05640454947398872     (/mV) 
  bn = -2.878962084018433     (1) 
  vhn = -81.91485211966854     (mV) 
  An = 0.3991683192025798     (/ms) 
  b1n = 0.038359520638365445     (/mV) 
  c1n = 0.0007470300267753836     (/mV2) 
  d1n = 1.2940312722420579e-05     (/mV3) 
  b2n = 0.03343351304321064     (/mV) 
  c2n = -0.0001494471549557467     (/mV2) 
  d2n = 2.508970605067338e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
  hInf 
  hTau 
  nInf 
  nTau 
}

STATE
{
  m
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*m*h
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  m = mInf 
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}