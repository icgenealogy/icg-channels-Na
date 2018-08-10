NEURON
{
  SUFFIX hh_wbm 
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

  am = 0.10469613962632877     (/mV) 
  bm = -3.6197993088871     (1) 
  vhm = -83.25000668373549     (mV) 
  Am = 0.026431768170894435     (/ms) 
  b1m = 0.08334309815112445     (/mV) 
  c1m = -0.0024102101428103494     (/mV2) 
  d1m = 1.3829257963919366e-05     (/mV3) 
  b2m = -6.864683550404374     (/mV) 
  c2m = 0.09130621537128153     (/mV2) 
  d2m = 0.000266918296118323     (/mV3) 

  ah = -0.14147535526867638     (/mV) 
  bh = 7.803714886955341     (1) 
  vhh = -59.842408425628626     (mV) 
  Ah = 0.12451629892151012     (/ms) 
  b1h = 0.08519686728812158     (/mV) 
  c1h = 0.0014408874589366281     (/mV2) 
  d1h = 1.301311301800623e-05     (/mV3) 
  b2h = 0.08647039911935184     (/mV) 
  c2h = -0.000895569613756051     (/mV2) 
  d2h = 2.8783176690042124e-06     (/mV3) 

  an = 0.057201680202445704     (/mV) 
  bn = -1.7152280884837368     (1) 
  vhn = -64.32199527937287     (mV) 
  An = 0.08083819385982961     (/ms) 
  b1n = 0.0003716891896213954     (/mV) 
  c1n = -0.0005415159573395605     (/mV2) 
  d1n = 3.7616971057495643e-06     (/mV3) 
  b2n = -0.00037238420122667895     (/mV) 
  c2n = -0.0005552569937489587     (/mV2) 
  d2n = 5.745762604832619e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
  hInf 
  hTau 
  nInf 
  nTau 
}

STATE
{
  m
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*m*h*n*n*n*n
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  m = mInf 
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}