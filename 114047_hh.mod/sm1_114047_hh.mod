NEURON
{
  SUFFIX HH 
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

  am = 0.10469647551301056     (/mV) 
  bm = -3.6198238021454805     (1) 
  vhm = -45.882586713780164     (mV) 
  Am = 0.3501206206037879     (/ms) 
  b1m = 0.021459093380793777     (/mV) 
  c1m = 0.0009791428224917484     (/mV2) 
  d1m = 2.8505979914657246e-05     (/mV3) 
  b2m = 3.9974674331226807     (/mV) 
  c2m = -0.0011284418976938238     (/mV2) 
  d2m = -7.801623737397019e-06     (/mV3) 

  ah = -0.14147587381380916     (/mV) 
  bh = 7.803745989974017     (1) 
  vhh = -45.04496626982896     (mV) 
  Ah = 0.37451091948443427     (/ms) 
  b1h = -0.11295845709698686     (/mV) 
  c1h = 0.0018065455525936248     (/mV2) 
  d1h = -7.5894696915615636e-06     (/mV3) 
  b2h = -0.002116107907523115     (/mV) 
  c2h = 0.0004894122958796934     (/mV2) 
  d2h = -2.744713300864651e-06     (/mV3) 

  an = 0.05720159481397284     (/mV) 
  bn = -1.715229123373379     (1) 
  vhn = -52.15867234439556     (mV) 
  An = 0.3987513599935466     (/ms) 
  b1n = -0.04119791493301328     (/mV) 
  c1n = 0.0002637015363414731     (/mV2) 
  d1n = -7.031312958568088e-07     (/mV3) 
  b2n = -0.03537453709688636     (/mV) 
  c2n = -0.0004997733369558299     (/mV2) 
  d2n = -4.308912534152548e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
  hInf 
  hTau 
  nInf 
  nTau 
}

STATE
{
  m
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*m*h*n*n*n*n
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  m = mInf 
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}