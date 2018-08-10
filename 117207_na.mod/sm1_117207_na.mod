NEURON
{
  SUFFIX na 
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

  ah = -0.16646448592306828     (/mV) 
  bh = 10.638092317238206     (1) 
  vhh = -69.70920562935412     (mV) 
  Ah = 2.17980113486916     (/ms) 
  b1h = -0.05258048015802101     (/mV) 
  c1h = 0.00036812127992216055     (/mV2) 
  d1h = -9.910539067292469e-07     (/mV3) 
  b2h = -0.10679766119241427     (/mV) 
  c2h = -0.0025430820383918733     (/mV2) 
  d2h = -3.082672214459087e-05     (/mV3) 

  am = 0.11105666566791125     (/mV) 
  bm = -4.270554876938829     (1) 
  vhm = -45.94568255222068     (mV) 
  Am = 0.23406055292377392     (/ms) 
  b1m = -0.041725556219116465     (/mV) 
  c1m = 0.00023266386832819475     (/mV2) 
  d1m = -4.394037448025439e-07     (/mV3) 
  b2m = -0.06469408128758541     (/mV) 
  c2m = -0.0010169499671365643     (/mV2) 
  d2m = -7.518718858915409e-06     (/mV3) 
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