NEURON
{
  SUFFIX Na 
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

  ah = -0.2474858777672662     (/mV) 
  bh = 9.483510177403815     (1) 
  vhh = -38.93505682688688     (mV) 
  Ah = 6.59410786547659     (/ms) 
  b1h = 0.06518504259725774     (/mV) 
  c1h = 0.0004131673084698631     (/mV2) 
  d1h = 4.695204183102238e-06     (/mV3) 
  b2h = 0.1988858743357794     (/mV) 
  c2h = -0.0027861871652129586     (/mV2) 
  d2h = 1.1485122196951142e-05     (/mV3) 

  am = 0.16139086559848378     (/mV) 
  bm = -5.655856641699078     (1) 
  vhm = -63.369797253162865     (mV) 
  Am = 0.11479681372127422     (/ms) 
  b1m = 0.0074782315176142485     (/mV) 
  c1m = -0.000496525515221604     (/mV2) 
  d1m = 2.757862564199573e-06     (/mV3) 
  b2m = -0.030614813169333167     (/mV) 
  c2m = -0.0008452728700796421     (/mV2) 
  d2m = 7.774039464724454e-06     (/mV3) 
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
  g = gbar*h*m*m
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