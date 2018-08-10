NEURON
{
  SUFFIX nap 
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

  am = 0.14999977976161896     (/mV) 
  bm = -8.399986801331593     (1) 
  vhm = -67.92745333338941     (mV) 
  Am = 10.26312539790677     (/ms) 
  b1m = -0.08823838267646326     (/mV) 
  c1m = 0.0009507592122651888     (/mV2) 
  d1m = -3.0863924792356612e-06     (/mV3) 
  b2m = -0.11012711070490049     (/mV) 
  c2m = -0.004000630170082858     (/mV2) 
  d2m = -6.499803873378804e-05     (/mV3) 
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
  g = gbar*m
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