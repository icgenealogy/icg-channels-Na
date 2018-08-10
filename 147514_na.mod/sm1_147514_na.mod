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

  ah = -0.16110669119942875     (/mV) 
  bh = 10.471940936899108     (1) 
  vhh = -66.93295278452013     (mV) 
  Ah = 19.214294548564535     (/ms) 
  b1h = -0.14048848410723316     (/mV) 
  c1h = 0.0014587842009776207     (/mV2) 
  d1h = -5.162286826814463e-06     (/mV3) 
  b2h = -0.14977270696326256     (/mV) 
  c2h = -0.00314913052480042     (/mV2) 
  d2h = -3.169011568855756e-05     (/mV3) 

  am = 0.11105651528115676     (/mV) 
  bm = -4.270554474418887     (1) 
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