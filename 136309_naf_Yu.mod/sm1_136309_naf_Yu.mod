NEURON
{
  SUFFIX nafYu 
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

  ah = -0.16127570922753418     (/mV) 
  bh = 8.870203340316712     (1) 
  vhh = -56.90188000288067     (mV) 
  Ah = 61.67335635997926     (/ms) 
  b1h = -0.1398558993194987     (/mV) 
  c1h = 0.0014038124256990395     (/mV2) 
  d1h = -4.63058449264023e-06     (/mV3) 
  b2h = -0.1473121864548035     (/mV) 
  c2h = -0.002867113637448324     (/mV2) 
  d2h = -2.4732409392967656e-05     (/mV3) 

  am = 0.11110884856309199     (/mV) 
  bm = -3.161454670401416     (1) 
  vhm = -26.846462216909877     (mV) 
  Am = 0.7511950455933766     (/ms) 
  b1m = 0.047693629644853595     (/mV) 
  c1m = 0.00039754535399106516     (/mV2) 
  d1m = 1.3676256020535733e-06     (/mV3) 
  b2m = 0.06077219727785108     (/mV) 
  c2m = -0.0005617206797559843     (/mV2) 
  d2m = 1.6837708027355717e-06     (/mV3) 
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