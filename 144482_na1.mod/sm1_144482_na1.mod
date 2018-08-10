NEURON
{
  SUFFIX na1 
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

  ah = -0.16112798854822644     (/mV) 
  bh = 9.909322043093155     (1) 
  vhh = -60.0     (mV) 
  Ah = 33.20000029854384     (/ms) 
  b1h = 13.861277367065137     (/mV) 
  c1h = -0.0016103186918793976     (/mV2) 
  d1h = 8.603952643820786e-06     (/mV3) 
  b2h = 13.859897336275774     (/mV) 
  c2h = 0.0013857875261616287     (/mV2) 
  d2h = -8.614875439059956e-06     (/mV3) 

  am = 0.11105507327671356     (/mV) 
  bm = -5.54756013106506     (1) 
  vhm = -59.35663637674245     (mV) 
  Am = 0.2232065577364239     (/ms) 
  b1m = -0.03851883295080149     (/mV) 
  c1m = 0.00018964937781525997     (/mV2) 
  d1m = -2.744958983974093e-07     (/mV3) 
  b2m = -0.07503011218723421     (/mV) 
  c2m = -0.001966074531934085     (/mV2) 
  d2m = -2.4777414465655495e-05     (/mV3) 
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