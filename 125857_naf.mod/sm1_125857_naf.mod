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

  ah = -0.24748682584546797     (/mV) 
  bh = 9.483605604842468     (1) 
  vhh = -38.982323184308946     (mV) 
  Ah = 10.369209909587559     (/ms) 
  b1h = -0.19888519422400902     (/mV) 
  c1h = 0.002783242909095485     (/mV2) 
  d1h = -1.1463340079946173e-05     (/mV3) 
  b2h = -0.0662910202143604     (/mV) 
  c2h = -0.0004598991997644932     (/mV2) 
  d2h = -5.282456317327823e-06     (/mV3) 
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