NEURON
{
  SUFFIX NaPyr 
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

  ah = -0.24749014437157882     (/mV) 
  bh = 11.21635023717488     (1) 
  vhh = -46.71468445546777     (mV) 
  Ah = 2.175949956577157     (/ms) 
  b1h = 0.07827221235606847     (/mV) 
  c1h = 0.0009802143157796345     (/mV2) 
  d1h = 1.1929034614265825e-05     (/mV3) 
  b2h = 0.19029410519162465     (/mV) 
  c2h = -0.002533645147696505     (/mV2) 
  d2h = 9.943064588503892e-06     (/mV3) 

  am = 0.1349135110247322     (/mV) 
  bm = -5.5348880555579285     (1) 
  vhm = 1.3632073322221612     (mV) 
  Am = 0.08014900478009696     (/ms) 
  b1m = 0.40134155523369824     (/mV) 
  c1m = 0.017880459664934324     (/mV2) 
  d1m = 0.00013901037990698406     (/mV3) 
  b2m = -0.003038945147920482     (/mV) 
  c2m = 0.0006696224268384745     (/mV2) 
  d2m = 5.96417149913295e-06     (/mV3) 
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