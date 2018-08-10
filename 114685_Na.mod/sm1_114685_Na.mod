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

  ah = -0.24748553074951626     (/mV) 
  bh = 9.483498999792651     (1) 
  vhh = -38.982323184308946     (mV) 
  Ah = 10.369209909587559     (/ms) 
  b1h = -0.19888519422400902     (/mV) 
  c1h = 0.002783242909095485     (/mV2) 
  d1h = -1.1463340079946173e-05     (/mV3) 
  b2h = -0.0662910202143604     (/mV) 
  c2h = -0.0004598991997644932     (/mV2) 
  d2h = -5.282456317327823e-06     (/mV3) 

  am = 0.13491522917982915     (/mV) 
  bm = -4.577042238553722     (1) 
  vhm = -55.16388616421631     (mV) 
  Am = 0.2004964610828142     (/ms) 
  b1m = -0.018021031880282777     (/mV) 
  c1m = -2.344299775078839e-05     (/mV2) 
  d1m = 3.91028635946108e-07     (/mV3) 
  b2m = -0.055869313975043466     (/mV) 
  c2m = -0.0009896312803736762     (/mV2) 
  d2m = -8.418419243792355e-06     (/mV3) 
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