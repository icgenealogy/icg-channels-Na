NEURON
{
  SUFFIX INap 
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

  ah = -0.12179938835028653     (/mV) 
  bh = 5.163958507336687     (1) 
  vhh = -45.94417014027623     (mV) 
  Ah = 4452.364992865754     (/ms) 
  b1h = 0.05712300779888295     (/mV) 
  c1h = 0.0008331342046289929     (/mV2) 
  d1h = 8.713455303463335e-06     (/mV3) 
  b2h = 0.0972788942119254     (/mV) 
  c2h = -0.0010727002185923946     (/mV2) 
  d2h = 3.7072873115096456e-06     (/mV3) 
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