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

  ah = -0.16110663687263654     (/mV) 
  bh = 8.860868737941011     (1) 
  vhh = -56.67113880665171     (mV) 
  Ah = 19.204871711417105     (/ms) 
  b1h = 0.14387703348475298     (/mV) 
  c1h = 0.0026910868730918567     (/mV2) 
  d1h = 2.2293860985077604e-05     (/mV3) 
  b2h = 0.14323290629690433     (/mV) 
  c2h = -0.0015409041303059722     (/mV2) 
  d2h = 5.686562359496021e-06     (/mV3) 

  am = 0.11105667365591557     (/mV) 
  bm = -3.159990502189752     (1) 
  vhm = -32.510796040273696     (mV) 
  Am = 0.24045947108446009     (/ms) 
  b1m = 0.0583938010876362     (/mV) 
  c1m = 0.0007129861490507317     (/mV2) 
  d1m = 3.713326850096412e-06     (/mV3) 
  b2m = 0.04813429209413229     (/mV) 
  c2m = -0.0003450886393518797     (/mV2) 
  d2m = 9.649154707784447e-07     (/mV3) 
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