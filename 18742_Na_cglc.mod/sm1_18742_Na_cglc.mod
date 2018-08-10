NEURON
{
  SUFFIX Na_cglc 
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

  ah = -0.21909779665514817     (/mV) 
  bh = 14.277091959674127     (1) 
  vhh = -63.71779823101979     (mV) 
  Ah = 209.0875329555687     (/ms) 
  b1h = 0.035295742433816     (/mV) 
  c1h = -0.0002331079478273922     (/mV2) 
  d1h = -2.5272716214430026e-07     (/mV3) 
  b2h = 0.20051038487141742     (/mV) 
  c2h = -0.0023626914841156223     (/mV2) 
  d2h = 8.217197900016886e-06     (/mV3) 

  am = 0.3390223821330544     (/mV) 
  bm = -13.56484613833178     (1) 
  vhm = -42.689887750972666     (mV) 
  Am = 7.420616773507917     (/ms) 
  b1m = -0.20145177012304147     (/mV) 
  c1m = 0.004259450743734857     (/mV2) 
  d1m = -2.958028759881412e-05     (/mV3) 
  b2m = -0.0859511223184971     (/mV) 
  c2m = -0.0022399790010931833     (/mV2) 
  d2m = -1.89157558526445e-05     (/mV3) 
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