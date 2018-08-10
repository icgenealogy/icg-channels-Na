NEURON
{
  SUFFIX NaxSH0_ChannelML 
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

  ah = -0.250001445559744     (/mV) 
  bh = 12.500082321346913     (1) 
  vhh = -48.583546034528744     (mV) 
  Ah = 15.022427173113703     (/ms) 
  b1h = 0.22572324233487864     (/mV) 
  c1h = 0.007147930818596841     (/mV2) 
  d1h = 7.85755973473173e-05     (/mV3) 
  b2h = 0.1985474798141924     (/mV) 
  c2h = -0.0027622906386054933     (/mV2) 
  d2h = 1.1032566521046917e-05     (/mV3) 

  am = 0.13888871974249817     (/mV) 
  bm = -5.337848831360608     (1) 
  vhm = -32.32194667743715     (mV) 
  Am = 0.239454431655582     (/ms) 
  b1m = 0.026899102240411036     (/mV) 
  c1m = -6.368411917764769e-05     (/mV2) 
  d1m = -1.925709597983284e-06     (/mV3) 
  b2m = 0.08162065171818465     (/mV) 
  c2m = -0.0008732029117506689     (/mV2) 
  d2m = 2.6256614816233587e-06     (/mV3) 
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
  g = gbar*h*m
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