NEURON
{
  SUFFIX nax 
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

  ah = -0.250000018841268     (/mV) 
  bh = 10.499994794247062     (1) 
  vhh = -40.51614655109287     (mV) 
  Ah = 15.654703934548925     (/ms) 
  b1h = -0.21383820569661988     (/mV) 
  c1h = 0.0031659712039694616     (/mV2) 
  d1h = -1.3454159447391807e-05     (/mV3) 
  b2h = -0.22928304449863718     (/mV) 
  c2h = -0.006849758435220738     (/mV2) 
  d2h = -6.887186710570096e-05     (/mV3) 

  am = 0.13888866959775245     (/mV) 
  bm = -4.226725438567228     (1) 
  vhm = -25.36767850366625     (mV) 
  Am = 0.24775962204894886     (/ms) 
  b1m = -0.07951075991333235     (/mV) 
  c1m = 0.0008218661080799935     (/mV2) 
  d1m = -2.286528072711291e-06     (/mV3) 
  b2m = -0.029154196575833912     (/mV) 
  c2m = 7.611234051296298e-05     (/mV2) 
  d2m = 2.4293431802569086e-06     (/mV3) 
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