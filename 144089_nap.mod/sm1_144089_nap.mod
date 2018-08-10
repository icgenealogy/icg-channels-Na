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

  ah = -0.25493047363855037     (/mV) 
  bh = 14.369960628083948     (1) 
  vhh = -56.364233394890796     (mV) 
  Ah = 1242.7644291470524     (/ms) 
  b1h = 0.24847806065277436     (/mV) 
  c1h = 2.5768883683825806e-06     (/mV2) 
  d1h = 5.476613895684293e-08     (/mV3) 
  b2h = 0.006467313397671033     (/mV) 
  c2h = -7.479625704817409e-07     (/mV2) 
  d2h = -2.8664152834323282e-09     (/mV3) 

  am = 0.09667127622806333     (/mV) 
  bm = -0.20064974238378133     (1) 
  vhm = 5.564668624952094     (mV) 
  Am = 0.2350652369655115     (/ms) 
  b1m = -0.02703762271942078     (/mV) 
  c1m = -8.880238646157239e-05     (/mV2) 
  d1m = 2.5402425191110227e-06     (/mV3) 
  b2m = -0.018963339231824233     (/mV) 
  c2m = 0.00013027974624251588     (/mV2) 
  d2m = 1.946959415994304e-06     (/mV3) 
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