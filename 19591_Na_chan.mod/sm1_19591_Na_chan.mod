NEURON
{
  SUFFIX naf_chan 
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

  ah = -0.17754873542508964     (/mV) 
  bh = 8.877439660174156     (1) 
  vhh = -45.402304936105     (mV) 
  Ah = 7.80860253348513     (/ms) 
  b1h = 0.05942184167348824     (/mV) 
  c1h = -0.00040875894806144007     (/mV2) 
  d1h = -2.190311722204831e-07     (/mV3) 
  b2h = 0.13126332320978984     (/mV) 
  c2h = -0.0014594021544309021     (/mV2) 
  d2h = 5.016402182433133e-06     (/mV3) 

  am = 0.14700020877586478     (/mV) 
  bm = -5.733019136882624     (1) 
  vhm = -37.568678817287264     (mV) 
  Am = 0.68185590069161     (/ms) 
  b1m = -0.10955994166927047     (/mV) 
  c1m = 0.0013752595590528554     (/mV2) 
  d1m = -5.204896561520481e-06     (/mV3) 
  b2m = -0.07257325491454501     (/mV) 
  c2m = -0.00023086499107008435     (/mV2) 
  d2m = 4.7528499990468895e-06     (/mV3) 
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