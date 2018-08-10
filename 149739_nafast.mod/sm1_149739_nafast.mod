NEURON
{
  SUFFIX nafast 
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

  ah = -0.2474993226546894     (/mV) 
  bh = 8.98914503942899     (1) 
  vhh = -36.79273853683692     (mV) 
  Ah = 10.250064854515056     (/ms) 
  b1h = 0.06339725804068175     (/mV) 
  c1h = 0.0003440747539987186     (/mV2) 
  d1h = 3.926689783721854e-06     (/mV3) 
  b2h = 0.2009163161850577     (/mV) 
  c2h = -0.0028527003530421794     (/mV2) 
  d2h = 1.1924785597701593e-05     (/mV3) 

  am = 0.13491572757365664     (/mV) 
  bm = -4.3207612360532694     (1) 
  vhm = -43.84257098088465     (mV) 
  Am = 0.23688187894101354     (/ms) 
  b1m = 0.053285372927072014     (/mV) 
  c1m = 0.0005480389691148648     (/mV2) 
  d1m = 1.664167768112946e-06     (/mV3) 
  b2m = 0.03246691416268017     (/mV) 
  c2m = -0.00015208663406051913     (/mV2) 
  d2m = 1.9100247172672962e-07     (/mV3) 
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