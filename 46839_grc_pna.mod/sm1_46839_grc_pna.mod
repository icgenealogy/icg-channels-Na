NEURON
{
  SUFFIX GrC_pNa 
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

  am = 0.19999956398885527     (/mV) 
  bm = -8.399973771255118     (1) 
  vhm = -48.13795606419214     (mV) 
  Am = 6.065163211980459     (/ms) 
  b1m = 0.11015664178807998     (/mV) 
  c1m = 0.0023426628458350918     (/mV2) 
  d1m = 2.0497336132892964e-05     (/mV3) 
  b2m = 0.07607213611390537     (/mV) 
  c2m = -0.0006894235557761777     (/mV2) 
  d2m = 2.3070536260501567e-06     (/mV3) 
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