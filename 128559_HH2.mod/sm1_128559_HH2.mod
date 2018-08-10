NEURON
{
  SUFFIX HH2 
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

  am = 0.13491526522800965     (/mV) 
  bm = -3.9159722099030563     (1) 
  vhm = -21.86612218125213     (mV) 
  Am = 0.20851430883133062     (/ms) 
  b1m = -0.04810030201469161     (/mV) 
  c1m = 0.0004427318448750986     (/mV2) 
  d1m = -1.4378231017334796e-06     (/mV3) 
  b2m = -0.025384453993543818     (/mV) 
  c2m = 6.37288706343224e-05     (/mV2) 
  d2m = 1.5780729122126667e-06     (/mV3) 

  ah = -0.24748823231570757     (/mV) 
  bh = 9.483622951389666     (1) 
  vhh = -38.97349331045796     (mV) 
  Ah = 9.28370627418893     (/ms) 
  b1h = -0.19879898007618788     (/mV) 
  c1h = 0.002786374346874089     (/mV2) 
  d1h = -1.1478538705172057e-05     (/mV3) 
  b2h = -0.0661110386532119     (/mV) 
  c2h = -0.0004537225633019906     (/mV2) 
  d2h = -5.234159146221094e-06     (/mV3) 

  an = 0.09068708828489183     (/mV) 
  bn = -2.6115064884088453     (1) 
  vhn = -51.29893418383732     (mV) 
  An = 2.886790249536415     (/ms) 
  b1n = -0.04191407570642935     (/mV) 
  c1n = 0.0002739553551948287     (/mV2) 
  d1n = -7.539535344990324e-07     (/mV3) 
  b2n = -0.06990263673383365     (/mV) 
  c2n = -0.0011921993227104525     (/mV2) 
  d2n = -1.096148589889475e-05     (/mV3) 
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
  g = gbar*m*m*m*h*n*n*n*n
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