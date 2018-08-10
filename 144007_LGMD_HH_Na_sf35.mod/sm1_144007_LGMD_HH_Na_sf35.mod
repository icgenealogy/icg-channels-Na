NEURON
{
  SUFFIX HH_Na35 
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

  ah = -0.1931401240090477     (/mV) 
  bh = 9.647210022895457     (1) 
  vhh = -43.01105293019582     (mV) 
  Ah = 1.0502849171535487     (/ms) 
  b1h = 0.03220852071409768     (/mV) 
  c1h = -0.0011046809754036973     (/mV2) 
  d1h = 5.680802530077598e-06     (/mV3) 
  b2h = 0.12434225493897508     (/mV) 
  c2h = -0.001862215525839453     (/mV2) 
  d2h = 7.729387414526134e-06     (/mV3) 

  am = 0.13703797650788418     (/mV) 
  bm = -5.4737472242846525     (1) 
  vhm = -79.31117659425789     (mV) 
  Am = 0.16486289328079037     (/ms) 
  b1m = -0.05243018998388957     (/mV) 
  c1m = 0.00026488840898675366     (/mV2) 
  d1m = 4.3999478659025126e-07     (/mV3) 
  b2m = -0.4918072725889055     (/mV) 
  c2m = 0.005698297941223684     (/mV2) 
  d2m = -1.630639148981715e-05     (/mV3) 
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