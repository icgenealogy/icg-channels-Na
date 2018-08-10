NEURON
{
  SUFFIX nap 
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

  ah = -0.2243564372448453     (/mV) 
  bh = 11.217596676854619     (1) 
  vhh = -123.29000819774934     (mV) 
  Ah = 7.969224082300975     (/ms) 
  b1h = 0.0027470868220483046     (/mV) 
  c1h = 2.057103346033686e-07     (/mV2) 
  d1h = 1.265817352747736e-07     (/mV3) 
  b2h = 0.0018551552585242886     (/mV) 
  c2h = 1.4033385901682472e-05     (/mV2) 
  d2h = -4.450740106340417e-08     (/mV3) 

  am = 0.16082452489809962     (/mV) 
  bm = -6.432963713937008     (1) 
  vhm = -49.99159813305387     (mV) 
  Am = 0.22409885232479146     (/ms) 
  b1m = -0.003362928345686276     (/mV) 
  c1m = 1.184380171927439e-05     (/mV2) 
  d1m = 1.6087391905163757e-06     (/mV3) 
  b2m = 0.0018793568050661713     (/mV) 
  c2m = 7.703599015479969e-05     (/mV2) 
  d2m = -3.873402945727699e-07     (/mV3) 
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
  g = gbar*h*h
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