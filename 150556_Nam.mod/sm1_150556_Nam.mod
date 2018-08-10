NEURON
{
  SUFFIX Nam 
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

  ah = -0.14147621654720113     (/mV) 
  bh = 6.813424600578755     (1) 
  vhh = -50.750601447829204     (mV) 
  Ah = 3.43156136883892     (/ms) 
  b1h = 0.07546604319544357     (/mV) 
  c1h = 0.0010750289944805112     (/mV2) 
  d1h = 1.3214292717073153e-05     (/mV3) 
  b2h = 0.09916394735275522     (/mV) 
  c2h = -0.0010785092032073076     (/mV2) 
  d2h = 3.679903935795175e-06     (/mV3) 

  am = 0.10469692656975246     (/mV) 
  bm = -2.88694642889346     (1) 
  vhm = -54.06363488034461     (mV) 
  Am = 3.426020965480995     (/ms) 
  b1m = 0.11439924211434824     (/mV) 
  c1m = 0.004998094774994443     (/mV2) 
  d1m = 0.0001230474722565586     (/mV3) 
  b2m = 0.08802309513439316     (/mV) 
  c2m = -0.0008390265073860967     (/mV2) 
  d2m = 2.519701448641648e-06     (/mV3) 
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
  g = gbar*h*m
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