NEURON
{
  SUFFIX na16 
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

  ah = -0.1611083109720216     (/mV) 
  bh = 10.794247932290169     (1) 
  vhh = -61.82368037141131     (mV) 
  Ah = 19.209935991416472     (/ms) 
  b1h = -0.14169417483125601     (/mV) 
  c1h = 0.001495622323227692     (/mV2) 
  d1h = -5.3982314794466665e-06     (/mV3) 
  b2h = -0.1471449286362768     (/mV) 
  c2h = -0.002915661716169059     (/mV2) 
  d2h = -2.644734492321099e-05     (/mV3) 

  am = 0.16647143645745366     (/mV) 
  bm = -6.709235504827171     (1) 
  vhm = -46.064781439286016     (mV) 
  Am = 0.3435596479514336     (/ms) 
  b1m = -0.06027300346545766     (/mV) 
  c1m = 0.0004582254358012716     (/mV2) 
  d1m = -1.24993643629886e-06     (/mV3) 
  b2m = -0.09084524781128508     (/mV) 
  c2m = -0.0015690535596363105     (/mV2) 
  d2m = -1.0444261572908726e-05     (/mV3) 
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