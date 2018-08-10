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

  ah = -0.22450865113924384     (/mV) 
  bh = 9.429423240990872     (1) 
  vhh = -46.034726130841634     (mV) 
  Ah = 4.670030529358142     (/ms) 
  b1h = 0.20327677031755614     (/mV) 
  c1h = 0.005796404428090007     (/mV2) 
  d1h = 5.2076900542167625e-05     (/mV3) 
  b2h = 0.10001100990841624     (/mV) 
  c2h = -0.0012570814717678144     (/mV2) 
  d2h = 4.729617922827223e-06     (/mV3) 

  am = 0.16088388107019175     (/mV) 
  bm = -6.114002641161545     (1) 
  vhm = -18.926177238645483     (mV) 
  Am = 0.12010780846042647     (/ms) 
  b1m = -0.001609071529615302     (/mV) 
  c1m = -0.00013147846399071648     (/mV2) 
  d1m = 8.641656961815854e-07     (/mV3) 
  b2m = 0.0630366127846457     (/mV) 
  c2m = -0.0018357846832602418     (/mV2) 
  d2m = 1.0286317099891602e-05     (/mV3) 
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
  g = gbar*h*h
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