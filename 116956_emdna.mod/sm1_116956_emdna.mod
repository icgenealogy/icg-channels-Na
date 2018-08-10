NEURON
{
  SUFFIX emdna 
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

  ah = -0.124999719919679     (/mV) 
  bh = 1.3749900921747336     (1) 
  vhh = -9.518513555892728     (mV) 
  Ah = 14.288447071622961     (/ms) 
  b1h = -0.08308845604405468     (/mV) 
  c1h = -1.9641859579576633e-05     (/mV2) 
  d1h = 4.00653623372086e-07     (/mV3) 
  b2h = -0.07092535178102714     (/mV) 
  c2h = 2.9062169278789222e-05     (/mV2) 
  d2h = 5.198924572374567e-07     (/mV3) 

  am = 0.16666650524244206     (/mV) 
  bm = -0.1666505288412264     (1) 
  vhm = -51.6765108563013     (mV) 
  Am = 3.9785532015234795     (/ms) 
  b1m = 0.09321055852123239     (/mV) 
  c1m = 8.033632129768526e-05     (/mV2) 
  d1m = 4.8185333835049e-07     (/mV3) 
  b2m = 0.03791116143217737     (/mV) 
  c2m = 3.8311894930754157e-05     (/mV2) 
  d2m = -2.723659659239226e-07     (/mV3) 
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