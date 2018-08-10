NEURON
{
  SUFFIX nakole 
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

  ah = -0.16110770419022563     (/mV) 
  bh = 8.860948199111242     (1) 
  vhh = -66.93295278452013     (mV) 
  Ah = 19.214294548564535     (/ms) 
  b1h = -0.14048848410723316     (/mV) 
  c1h = 0.0014587842009776207     (/mV2) 
  d1h = -5.162286826814463e-06     (/mV3) 
  b2h = -0.14977270696326256     (/mV) 
  c2h = -0.00314913052480042     (/mV2) 
  d2h = -3.169011568855756e-05     (/mV3) 

  am = 0.11105869910809782     (/mV) 
  bm = -3.4932363200970404     (1) 
  vhm = -40.480161764443764     (mV) 
  Am = 0.23597607091643272     (/ms) 
  b1m = -0.042796244236841016     (/mV) 
  c1m = 0.00025144444962337225     (/mV2) 
  d1m = -5.298162613710425e-07     (/mV3) 
  b2m = -0.06419477430664235     (/mV) 
  c2m = -0.0008889099506394048     (/mV2) 
  d2m = -5.196687886285116e-06     (/mV3) 
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