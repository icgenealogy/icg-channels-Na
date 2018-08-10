NEURON
{
  SUFFIX hhI 
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

  ah = -0.14147535022832108     (/mV) 
  bh = 7.80371458314383     (1) 
  vhh = -59.842408425628626     (mV) 
  Ah = 0.12451629892151012     (/ms) 
  b1h = 0.08519686728812158     (/mV) 
  c1h = 0.0014408874589366281     (/mV2) 
  d1h = 1.301311301800623e-05     (/mV3) 
  b2h = 0.08647039911935184     (/mV) 
  c2h = -0.000895569613756051     (/mV2) 
  d2h = 2.8783176690042124e-06     (/mV3) 

  an = 0.05720154328570853     (/mV) 
  bn = -1.71522199197761     (1) 
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
  hInf 
  hTau 
  nInf 
  nTau 
}

STATE
{
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*n*n*n*n
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}