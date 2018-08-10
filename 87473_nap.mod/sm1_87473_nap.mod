NEURON
{
  SUFFIX nap 
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

  ap = 0.1500000014379106     (/mV) 
  bp = -8.400001150422652     (1) 
  vhp = -73.56183466673515     (mV) 
  Ap = 10.020000007078417     (/ms) 
  b1p = -3.322440488619548e-09     (/mV) 
  c1p = 1.803021258416065e-09     (/mV2) 
  d1p = -5.310905572541754e-11     (/mV3) 
  b2p = -3.3883273721517777e-09     (/mV) 
  c2p = 1.8063835814949279e-09     (/mV2) 
  d2p = -5.3134583518005927e-11     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  pInf 
  pTau 
}

STATE
{
  p
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*p
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  p' = (pInf - p) / pTau 
}

INITIAL
{
  rates(v)
  p = pInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    pInf = 1/(1 + exp(-ap*v + bp)) 
    pTau = Ap / ( exp(-(b1p*(v-vhp) + c1p*(v-vhp)^2 + d1p*(v-vhp)^3)) + exp((b2p*(v-vhp) + c2p*(v-vhp)^2 + d2p*(v-vhp)^3)) ) 


  UNITSON
}