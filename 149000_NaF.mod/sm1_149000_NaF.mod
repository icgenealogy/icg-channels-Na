NEURON
{
  SUFFIX NaF 
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

  ah = -0.11195126879439073     (/mV) 
  bh = 8.64148802539199     (1) 
  vhh = -79.25703931074575     (mV) 
  Ah = 10.75749701214114     (/ms) 
  b1h = 0.060510999525982746     (/mV) 
  c1h = 0.001692444842009529     (/mV2) 
  d1h = 1.8045360460208446e-05     (/mV3) 
  b2h = 0.04673169984196739     (/mV) 
  c2h = 0.0001373947396536091     (/mV2) 
  d2h = -7.223843407670424e-07     (/mV3) 

  am = 0.15000099535521802     (/mV) 
  bm = -5.3595106101320065     (1) 
  vhm = -26.724276857026055     (mV) 
  Am = 0.4054320947104264     (/ms) 
  b1m = 0.0013343518679176907     (/mV) 
  c1m = -0.000755570156385375     (/mV2) 
  d1m = -3.454238423310397e-06     (/mV3) 
  b2m = 0.17298544541263056     (/mV) 
  c2m = -0.0025163685537788363     (/mV2) 
  d2m = -6.680018782707187e-05     (/mV3) 
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