NEURON
{
  SUFFIX napf_spinstell 
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

  am = 0.10017344991838191     (/mV) 
  bm = -3.558508436363909     (1) 
  vhm = -30.55179893806154     (mV) 
  Am = 0.2693685431147677     (/ms) 
  b1m = -0.10476719715925986     (/mV) 
  c1m = 0.001435480372812827     (/mV2) 
  d1m = -5.95202459899755e-06     (/mV3) 
  b2m = -0.16323350228321107     (/mV) 
  c2m = -0.0038822406531991     (/mV2) 
  d2m = -2.9150364163998165e-05     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*m
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}