NEURON
{
  SUFFIX hh2 
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

  ah = -0.24749376811385332     (/mV) 
  bh = 10.226329175300414     (1) 
  vhh = -42.173372065747934     (mV) 
  Ah = 9.405292087036221     (/ms) 
  b1h = 0.06940075318597035     (/mV) 
  c1h = 0.0005847610715628065     (/mV2) 
  d1h = 6.800648326486456e-06     (/mV3) 
  b2h = 0.19691864303689924     (/mV) 
  c2h = -0.0027073183640460036     (/mV2) 
  d2h = 1.0922774529196114e-05     (/mV3) 

  an = 0.09068018681921858     (/mV) 
  bn = -3.3367437295561415     (1) 
  vhn = -60.05724641942809     (mV) 
  An = 2.8448625738862363     (/ms) 
  b1n = 0.0711545854985429     (/mV) 
  c1n = 0.001353703424339468     (/mV2) 
  d1n = 1.4586364001697535e-05     (/mV3) 
  b2n = 0.04018794278995792     (/mV) 
  c2n = -0.00024655570589804865     (/mV2) 
  d2n = 6.373542640127566e-07     (/mV3) 

  am = 0.13491514523045064     (/mV) 
  bm = -4.995277507004642     (1) 
  vhm = -42.04952591849492     (mV) 
  Am = 0.2228202428937888     (/ms) 
  b1m = -0.03846787226732772     (/mV) 
  c1m = 0.00023272998052925107     (/mV2) 
  d1m = -4.207735166708235e-07     (/mV3) 
  b2m = -0.04060140157540383     (/mV) 
  c2m = -7.648309385394124e-05     (/mV2) 
  d2m = 3.1336141634278755e-06     (/mV3) 
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
  nInf 
  nTau 
  mInf 
  mTau 
}

STATE
{
  h
  n
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*n*n*n*n*m*m*m
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  n = nInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}