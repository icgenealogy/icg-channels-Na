NEURON
{
  SUFFIX nav1p8 
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

  ah = -0.10413261320467776     (/mV) 
  bh = 4.228702065371414     (1) 
  vhh = -37.23751953126361     (mV) 
  Ah = 263.3088631954705     (/ms) 
  b1h = 0.011183890769101704     (/mV) 
  c1h = -0.00019452005144610261     (/mV2) 
  d1h = -1.2340201045796645e-06     (/mV3) 
  b2h = 0.09088930603441421     (/mV) 
  c2h = -0.00013044313325846242     (/mV2) 
  d2h = -2.5697753161667954e-06     (/mV3) 

  am = 0.11055741262593959     (/mV) 
  bm = -2.0605242911287016     (1) 
  vhm = -26.477814863842333     (mV) 
  Am = 1.401038081248401     (/ms) 
  b1m = -0.06138071390693101     (/mV) 
  c1m = 0.0007249483540798967     (/mV2) 
  d1m = -2.732044835947465e-06     (/mV3) 
  b2m = -0.0715937544766443     (/mV) 
  c2m = -0.0009285849882232505     (/mV2) 
  d2m = -4.760008123342363e-06     (/mV3) 
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