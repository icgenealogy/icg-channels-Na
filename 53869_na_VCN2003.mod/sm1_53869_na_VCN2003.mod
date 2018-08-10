NEURON
{
  SUFFIX na_VCN2003 
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

  ah = -0.16666655836082878     (/mV) 
  bh = 10.833329618174313     (1) 
  vhh = -58.67743765657265     (mV) 
  Ah = 2.4695776703393153     (/ms) 
  b1h = 0.05927806377050525     (/mV) 
  c1h = 0.0009457540502650657     (/mV2) 
  d1h = 1.0750382498636077e-05     (/mV3) 
  b2h = 0.09745928879619518     (/mV) 
  c2h = -0.00097578532678134     (/mV2) 
  d2h = 3.076397238760479e-06     (/mV3) 

  am = 0.14285677935534838     (/mV) 
  bm = -5.428558277215727     (1) 
  vhm = -50.17200724898032     (mV) 
  Am = 0.13961873499575223     (/ms) 
  b1m = -0.025113159486536446     (/mV) 
  c1m = -0.00025106065406131455     (/mV2) 
  d1m = 2.4323832011467046e-06     (/mV3) 
  b2m = -0.04242150969355313     (/mV) 
  c2m = -5.455189800687058e-05     (/mV2) 
  d2m = 2.640359398796647e-06     (/mV3) 
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
  g = gbar*h*m*m
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