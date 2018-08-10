NEURON
{
  SUFFIX fna 
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

  ah = -0.1459111291475059     (/mV) 
  bh = 7.230420551740632     (1) 
  vhh = -46.16320459454851     (mV) 
  Ah = 9.280314795320958     (/ms) 
  b1h = 0.03516328295847089     (/mV) 
  c1h = -0.00019894948126867711     (/mV2) 
  d1h = -2.576142806652954e-07     (/mV3) 
  b2h = 0.12108111249928602     (/mV) 
  c2h = -0.0013294636755214688     (/mV2) 
  d2h = 4.569752482838888e-06     (/mV3) 

  am = 0.12492406314908665     (/mV) 
  bm = -3.8726501945967406     (1) 
  vhm = -46.15666850231317     (mV) 
  Am = 0.21207985447459324     (/ms) 
  b1m = 0.055494283792121936     (/mV) 
  c1m = 0.0006366182161894302     (/mV2) 
  d1m = 1.6951087910904022e-06     (/mV3) 
  b2m = 0.024862516404061175     (/mV) 
  c2m = -7.817066959664852e-05     (/mV2) 
  d2m = 6.825642217128833e-09     (/mV3) 
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