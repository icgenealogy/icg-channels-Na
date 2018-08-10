NEURON
{
  SUFFIX naxn 
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

  ah = -0.25000024578415325     (/mV) 
  bh = 12.500014163778847     (1) 
  vhh = -48.58580744292182     (mV) 
  Ah = 15.021436940061905     (/ms) 
  b1h = -0.19848284370367575     (/mV) 
  c1h = 0.002761663795063771     (/mV2) 
  d1h = -1.103186821249092e-05     (/mV3) 
  b2h = -0.22576001227361392     (/mV) 
  c2h = -0.007150300413447884     (/mV2) 
  d2h = -7.861107637880668e-05     (/mV3) 

  am = 0.13888876345689102     (/mV) 
  bm = -5.337842997485151     (1) 
  vhm = -33.22076321875721     (mV) 
  Am = 0.24453680855779153     (/ms) 
  b1m = 0.02872048633845459     (/mV) 
  c1m = -4.8858328561324556e-05     (/mV2) 
  d1m = -2.031721975744898e-06     (/mV3) 
  b2m = 0.07908667846228963     (/mV) 
  c2m = -0.0007873040389808289     (/mV2) 
  d2m = 2.2065926644326795e-06     (/mV3) 
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