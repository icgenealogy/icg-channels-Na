NEURON
{
  SUFFIX naf 
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

  ah = -0.14147536116666404     (/mV) 
  bh = 7.803715242145243     (1) 
  vhh = -59.02087970052934     (mV) 
  Ah = 3.470017381800921     (/ms) 
  b1h = -0.0952946367142484     (/mV) 
  c1h = 0.000988338842643844     (/mV2) 
  d1h = -3.2100820268035223e-06     (/mV3) 
  b2h = -0.08646919305677003     (/mV) 
  c2h = -0.0016798615376884519     (/mV2) 
  d2h = -2.3507338091177925e-05     (/mV3) 
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
}

STATE
{
  h
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
}

INITIAL
{
  rates(v)
  h = hInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 


  UNITSON
}