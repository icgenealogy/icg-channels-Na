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

  ah = -0.2670730190281382     (/mV) 
  bh = 10.197315116576378     (1) 
  vhh = -45.19226418714471     (mV) 
  Ah = 12.921509812968434     (/ms) 
  b1h = -0.12884369468945378     (/mV) 
  c1h = 0.001379004619609996     (/mV2) 
  d1h = -5.184714827238802e-06     (/mV3) 
  b2h = -0.04606318893287176     (/mV) 
  c2h = -0.0006149408271058172     (/mV2) 
  d2h = -4.189835202501841e-06     (/mV3) 

  am = 0.10731362895519561     (/mV) 
  bm = -2.5393209095438807     (1) 
  vhm = -30.768146936845717     (mV) 
  Am = 0.10276099205845286     (/ms) 
  b1m = -0.0020382272741375387     (/mV) 
  c1m = -0.00022695978048671572     (/mV2) 
  d1m = 1.445991575265036e-06     (/mV3) 
  b2m = -0.043709459628985786     (/mV) 
  c2m = -0.00025048225674579705     (/mV2) 
  d2m = 4.545693830138169e-06     (/mV3) 
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