NEURON
{
  SUFFIX nargc 
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

  ah = -0.14209086032953172     (/mV) 
  bh = 6.818399595988758     (1) 
  vhh = -49.91762006868555     (mV) 
  Ah = 6.546128085964046     (/ms) 
  b1h = 0.07081361876028652     (/mV) 
  c1h = 0.0008862215815560009     (/mV2) 
  d1h = 1.0758443771341796e-05     (/mV3) 
  b2h = 0.10175344568422655     (/mV) 
  c2h = -0.001099681072969431     (/mV2) 
  d2h = 3.722236828794977e-06     (/mV3) 

  am = 0.15877746943459756     (/mV) 
  bm = -4.835363633713776     (1) 
  vhm = -37.11577311289199     (mV) 
  Am = 0.43411158200713246     (/ms) 
  b1m = -0.07072341846062224     (/mV) 
  c1m = 0.000683170962068082     (/mV2) 
  d1m = -2.4918694288728364e-06     (/mV3) 
  b2m = -0.11316464129943841     (/mV) 
  c2m = -0.0019387888253864667     (/mV2) 
  d2m = -2.119956888286355e-05     (/mV3) 
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