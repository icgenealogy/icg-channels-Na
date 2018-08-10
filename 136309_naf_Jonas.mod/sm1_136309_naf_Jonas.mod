NEURON
{
  SUFFIX nafJonas 
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

  ah = -0.1281093732363397     (/mV) 
  bh = 9.5748547974122     (1) 
  vhh = -68.74031177270254     (mV) 
  Ah = 23.911423984655322     (/ms) 
  b1h = -0.0939765460885211     (/mV) 
  c1h = 0.00040068217037145313     (/mV2) 
  d1h = 1.9736094043841233e-07     (/mV3) 
  b2h = -0.026816409063559685     (/mV) 
  c2h = 0.00040911210566089843     (/mV2) 
  d2h = -4.188853543164363e-07     (/mV3) 

  am = 0.09312496605572201     (/mV) 
  bm = -4.00168088537632     (1) 
  vhm = 42.65868688437798     (mV) 
  Am = 0.9241217419206534     (/ms) 
  b1m = -1.249433887872804     (/mV) 
  c1m = -0.025363741707544116     (/mV2) 
  d1m = -0.00012802500245441735     (/mV3) 
  b2m = -0.2220054943367586     (/mV) 
  c2m = -0.004422573684507313     (/mV2) 
  d2m = -2.28401283604296e-05     (/mV3) 
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