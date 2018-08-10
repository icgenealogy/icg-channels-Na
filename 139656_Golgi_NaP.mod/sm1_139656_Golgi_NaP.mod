NEURON
{
  SUFFIX Golgi_NaP 
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

  am = 0.19999998714708458     (/mV) 
  bm = -8.599994617333158     (1) 
  vhm = -45.51923818135567     (mV) 
  Am = 0.6303343739315928     (/ms) 
  b1m = 0.10410501317059219     (/mV) 
  c1m = 0.0020745686264167854     (/mV2) 
  d1m = 1.6766597625551468e-05     (/mV3) 
  b2m = 0.07873074982066325     (/mV) 
  c2m = -0.0007879132304048668     (/mV2) 
  d2m = 2.542103678328071e-06     (/mV3) 
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
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}