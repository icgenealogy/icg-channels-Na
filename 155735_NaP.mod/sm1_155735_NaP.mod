NEURON
{
  SUFFIX NaP 
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

  am = 0.15989849753542346     (/mV) 
  bm = -7.177593054170288     (1) 
  vhm = -41.64967064031089     (mV) 
  Am = 0.27700050982841373     (/ms) 
  b1m = -0.10022496246866677     (/mV) 
  c1m = 0.0017318870499612364     (/mV2) 
  d1m = -1.508402205210748e-05     (/mV3) 
  b2m = -0.10895738605079291     (/mV) 
  c2m = -0.0022191778718422908     (/mV2) 
  d2m = -1.4575243423466424e-05     (/mV3) 
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