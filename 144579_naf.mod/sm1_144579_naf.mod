NEURON
{
  SUFFIX naf 
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

  ah = -0.09345771429840262     (/mV) 
  bh = 5.551402752750402     (1) 
  vhh = -25.453236181644463     (mV) 
  Ah = 1.1482175967605177     (/ms) 
  b1h = 0.0026052509864933     (/mV) 
  c1h = 0.00014367974258837511     (/mV2) 
  d1h = 1.18161502974378e-06     (/mV3) 
  b2h = 0.0644573356915108     (/mV) 
  c2h = -0.0006993244471533419     (/mV2) 
  d2h = 2.4931349758865164e-06     (/mV3) 

  am = 0.09999945658440339     (/mV) 
  bm = -3.45009272203922     (1) 
  vhm = -29.340728753887635     (mV) 
  Am = 0.27021781884301255     (/ms) 
  b1m = -0.10589059625219208     (/mV) 
  c1m = 0.0014665923573438685     (/mV2) 
  d1m = -6.1482308615545405e-06     (/mV3) 
  b2m = -0.15963561110762725     (/mV) 
  c2m = -0.00367516319767853     (/mV2) 
  d2m = -2.670392235092432e-05     (/mV3) 
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