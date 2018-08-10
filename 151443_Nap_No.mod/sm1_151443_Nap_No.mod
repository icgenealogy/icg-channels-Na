NEURON
{
  SUFFIX Nap_No 
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

  ap = 0.10164684435139372     (/mV) 
  bp = -5.741086555513508     (1) 
  vhp = -49.54225951622753     (mV) 
  Ap = 27.126493328591007     (/ms) 
  b1p = 0.010436833193929882     (/mV) 
  c1p = -6.224093978319443e-05     (/mV2) 
  d1p = -5.751059773475771e-07     (/mV3) 
  b2p = 0.08600853068046764     (/mV) 
  c2p = -0.0006107874149555077     (/mV2) 
  d2p = 1.6422277404591437e-06     (/mV3) 
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
  g = gbar*p*p*p
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