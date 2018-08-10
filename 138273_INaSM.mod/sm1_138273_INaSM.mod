NEURON
{
  SUFFIX INaSM 
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

  ah = -0.16128859788516725     (/mV) 
  bh = 10.483810710044308     (1) 
  vhh = -66.91295241475018     (mV) 
  Ah = 61.74642015149949     (/ms) 
  b1h = -0.1409599417947419     (/mV) 
  c1h = 0.001466422421736479     (/mV2) 
  d1h = -5.206522705589487e-06     (/mV3) 
  b2h = -0.14978235637868867     (/mV) 
  c2h = -0.003141691058317135     (/mV2) 
  d2h = -3.1627034405983603e-05     (/mV3) 

  am = 0.10636040164634604     (/mV) 
  bm = -4.372644955678366     (1) 
  vhm = -45.756166345694645     (mV) 
  Am = 0.6592939953185941     (/ms) 
  b1m = 0.05332400204984347     (/mV) 
  c1m = 0.0005592923072945726     (/mV2) 
  d1m = 2.3221217648841804e-06     (/mV3) 
  b2m = 0.04577870518631979     (/mV) 
  c2m = -0.00030210081880244124     (/mV2) 
  d2m = 7.940287877700027e-07     (/mV3) 
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