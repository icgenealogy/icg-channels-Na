NEURON
{
  SUFFIX naf2 
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

  ah = -0.1492475050926597     (/mV) 
  bh = 8.701154690771363     (1) 
  vhh = -48.64880354200927     (mV) 
  Ah = 2.0054516445282937     (/ms) 
  b1h = -0.06113328266795243     (/mV) 
  c1h = 0.0005688956009677608     (/mV2) 
  d1h = -1.7158293004390988e-06     (/mV3) 
  b2h = -0.027988358185658115     (/mV) 
  c2h = -0.0006600690834449608     (/mV2) 
  d2h = -5.205725054857754e-06     (/mV3) 

  am = 0.09999662402367916     (/mV) 
  bm = -3.799923108946324     (1) 
  vhm = -23.04576891218555     (mV) 
  Am = 0.2186617300056737     (/ms) 
  b1m = -0.14859158941137357     (/mV) 
  c1m = 0.0031634880694006487     (/mV2) 
  d1m = -1.6809276423736407e-05     (/mV3) 
  b2m = -0.02883114258311653     (/mV) 
  c2m = 0.0012782759777738696     (/mV2) 
  d2m = -7.964195760349671e-06     (/mV3) 
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