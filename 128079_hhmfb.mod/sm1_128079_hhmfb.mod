NEURON
{
  SUFFIX hhmfb 
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

  am = 0.09314983481016716     (/mV) 
  bm = -4.001348709275805     (1) 
  vhm = -24.870972633767177     (mV) 
  Am = 0.40524004812522235     (/ms) 
  b1m = 0.0031430003364563783     (/mV) 
  c1m = -0.000253368754268289     (/mV2) 
  d1m = 1.3207459860424574e-06     (/mV3) 
  b2m = 0.07462424284024048     (/mV) 
  c2m = 3.8727007326569354e-05     (/mV2) 
  d2m = -1.1642489689961596e-05     (/mV3) 

  ah = -0.12803706304853585     (/mV) 
  bh = 9.574697132102136     (1) 
  vhh = -68.74929735981553     (mV) 
  Ah = 23.89382161746821     (/ms) 
  b1h = -0.09401574166602678     (/mV) 
  c1h = 0.00040246775277203513     (/mV2) 
  d1h = 1.874859814697124e-07     (/mV3) 
  b2h = -0.026665846281499803     (/mV) 
  c2h = 0.0004111522458296057     (/mV2) 
  d2h = -4.3773584865741607e-07     (/mV3) 

  an = 0.05640443274313314     (/mV) 
  bn = -2.8789538938026533     (1) 
  vhn = -78.75989675347725     (mV) 
  An = 11.580293927518788     (/ms) 
  b1n = 0.038871524726985254     (/mV) 
  c1n = 0.0005648178175757243     (/mV2) 
  d1n = 6.075993390528273e-06     (/mV3) 
  b2n = 0.036592813882797906     (/mV) 
  c2n = -0.0001791165122602799     (/mV2) 
  d2n = 3.6493932826050375e-07     (/mV3) 
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