NEURON
{
  SUFFIX Nap 
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

  ah = -0.09964518934198312     (/mV) 
  bh = 4.866934780794049     (1) 
  vhh = -7.648353442027916     (mV) 
  Ah = 1276.4326223188025     (/ms) 
  b1h = -0.014606280129179555     (/mV) 
  c1h = 0.0003600529233861533     (/mV2) 
  d1h = -2.5739549530966196e-06     (/mV3) 
  b2h = 0.014594601651732571     (/mV) 
  c2h = 9.087432507501431e-05     (/mV2) 
  d2h = -1.0121242691439518e-06     (/mV3) 
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