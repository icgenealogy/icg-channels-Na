NEURON
{
  SUFFIX hh3 
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

  am = 0.33335587823103086     (/mV) 
  bm = -13.335675879046082     (1) 
  vhm = -219.9321323598076     (mV) 
  Am = 0.00017938862572965602     (/ms) 
  b1m = 0.11107262861596738     (/mV) 
  c1m = -0.0006858270702761768     (/mV2) 
  d1m = 1.3628702341933297e-06     (/mV3) 
  b2m = -0.1776837136689594     (/mV) 
  c2m = 0.0009653948124177064     (/mV2) 
  d2m = -1.4544559457615275e-06     (/mV3) 

  ah = -0.3334685379700192     (/mV) 
  bh = 15.006778011378302     (1) 
  vhh = -303.72889470996375     (mV) 
  Ah = 0.3274152235454768     (/ms) 
  b1h = 0.004526652172878978     (/mV) 
  c1h = -1.4610367430325444e-05     (/mV2) 
  d1h = 1.5164698776552223e-08     (/mV3) 
  b2h = -0.2294905565551826     (/mV) 
  c2h = 0.0014948372429448922     (/mV2) 
  d2h = -2.58431030337447e-06     (/mV3) 

  an = 0.3333558782435763     (/mV) 
  bn = -13.335675879574337     (1) 
  vhn = -45.505977729033404     (mV) 
  An = 4.005311148887427     (/ms) 
  b1n = -0.00023109748976549307     (/mV) 
  c1n = 3.489663111228414e-05     (/mV2) 
  d1n = -2.4027910641749896e-07     (/mV3) 
  b2n = -7.172389471874845e-05     (/mV) 
  c2n = 3.044172756831485e-05     (/mV2) 
  d2n = -2.1965358430605933e-07     (/mV3) 
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
  g = gbar*m*m*h*n*n
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