NEURON
{
  SUFFIX NaTs2_t 
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

  ah = -0.16666644481154583     (/mV) 
  bh = 9.999985670077569     (1) 
  vhh = -38.188961700579625     (mV) 
  Ah = 2.0650132234858765     (/ms) 
  b1h = -0.07452367932585287     (/mV) 
  c1h = 0.0013199208325084808     (/mV2) 
  d1h = -6.322989013993274e-06     (/mV3) 
  b2h = 0.006846541105555631     (/mV) 
  c2h = 0.00039084084166815026     (/mV2) 
  d2h = -2.3522368544644017e-06     (/mV3) 

  am = 0.16666647540131668     (/mV) 
  bm = -5.717043938699035     (1) 
  vhm = -37.42172233816433     (mV) 
  Am = 0.3817380340926428     (/ms) 
  b1m = 0.08284445756459356     (/mV) 
  c1m = 0.0013239772763240607     (/mV2) 
  d1m = 8.584983383992002e-06     (/mV3) 
  b2m = 0.06973705145796948     (/mV) 
  c2m = -0.000662081641244349     (/mV2) 
  d2m = 2.340172616660308e-06     (/mV3) 
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