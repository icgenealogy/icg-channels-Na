NEURON
{
  SUFFIX na 
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

  ah = -0.1611085311772695     (/mV) 
  bh = 10.794263854506594     (1) 
  vhh = -61.82368037141131     (mV) 
  Ah = 19.209935991416472     (/ms) 
  b1h = -0.14169417483125601     (/mV) 
  c1h = 0.001495622323227692     (/mV2) 
  d1h = -5.3982314794466665e-06     (/mV3) 
  b2h = -0.1471449286362768     (/mV) 
  c2h = -0.002915661716169059     (/mV2) 
  d2h = -2.644734492321099e-05     (/mV3) 

  am = 0.14273535732516596     (/mV) 
  bm = -5.807382291622256     (1) 
  vhm = -45.79142554012623     (mV) 
  Am = 0.2977622451254423     (/ms) 
  b1m = -0.05614518226414453     (/mV) 
  c1m = 0.00041742780829248687     (/mV2) 
  d1m = -1.1158730244024082e-06     (/mV3) 
  b2m = -0.07514624526427317     (/mV) 
  c2m = -0.0010534977443832914     (/mV2) 
  d2m = -5.231890566595067e-06     (/mV3) 
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