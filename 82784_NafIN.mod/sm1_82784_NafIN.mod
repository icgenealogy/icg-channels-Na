NEURON
{
  SUFFIX NafIN 
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

  ah = -0.1414763111515706     (/mV) 
  bh = 7.110533928058724     (1) 
  vhh = -40.28963001933216     (mV) 
  Ah = 3.877239802722704     (/ms) 
  b1h = -0.11404914432410178     (/mV) 
  c1h = 0.001744829125061949     (/mV2) 
  d1h = -7.153745700746437e-06     (/mV3) 
  b2h = -0.0030423635345090115     (/mV) 
  c2h = 0.00047566560620845876     (/mV2) 
  d2h = -2.530780239162846e-06     (/mV3) 

  am = 0.09892000743461596     (/mV) 
  bm = -2.88110859832115     (1) 
  vhm = -36.33889368205952     (mV) 
  Am = 0.25198361374518025     (/ms) 
  b1m = 0.045129932721113236     (/mV) 
  c1m = 0.0004925722980072289     (/mV2) 
  d1m = 2.1900098982108527e-06     (/mV3) 
  b2m = 0.028761083519054906     (/mV) 
  c2m = -0.00012290554822013986     (/mV2) 
  d2m = 8.816801050769568e-08     (/mV3) 
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