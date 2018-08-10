NEURON
{
  SUFFIX nattxs 
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

  ah = -0.12588383352225718     (/mV) 
  bh = 9.179141159751161     (1) 
  vhh = -71.77375850939856     (mV) 
  Ah = 151.1859101674789     (/ms) 
  b1h = -0.07909494894664883     (/mV) 
  c1h = -2.9378187503557378e-05     (/mV2) 
  d1h = 1.967090946017215e-06     (/mV3) 
  b2h = -0.04277933093478552     (/mV) 
  c2h = 0.00018428233874542466     (/mV2) 
  d2h = 1.685804460098336e-06     (/mV3) 

  am = 0.11808318177022932     (/mV) 
  bm = -2.546958957947128     (1) 
  vhm = -32.043654288560006     (mV) 
  Am = 0.6313180329566356     (/ms) 
  b1m = -0.07840715507673379     (/mV) 
  c1m = 0.000963351859433003     (/mV2) 
  d1m = -3.6920610947577844e-06     (/mV3) 
  b2m = -0.06737824620417097     (/mV) 
  c2m = -0.0011869026261867636     (/mV2) 
  d2m = -8.380116248897237e-06     (/mV3) 
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