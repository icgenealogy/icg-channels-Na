NEURON
{
  SUFFIX nat 
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

  ah = -0.16122399134013538     (/mV) 
  bh = 8.867514021813943     (1) 
  vhh = -66.92696780998283     (mV) 
  Ah = 19.226711224616743     (/ms) 
  b1h = 0.14976631476087107     (/mV) 
  c1h = 0.0031382479334505094     (/mV2) 
  d1h = 3.143649695024415e-05     (/mV3) 
  b2h = 0.1406178445187404     (/mV) 
  c2h = -0.0014617778842253252     (/mV2) 
  d2h = 5.178683707978736e-06     (/mV3) 

  am = 0.11109466911544436     (/mV) 
  bm = -3.494414162221607     (1) 
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