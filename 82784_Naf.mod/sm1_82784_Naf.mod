NEURON
{
  SUFFIX Naf 
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

  ah = -0.1414762511308444     (/mV) 
  bh = 5.695767112484065     (1) 
  vhh = -41.461713164019386     (mV) 
  Ah = 11.955601177493413     (/ms) 
  b1h = 0.06437014136701824     (/mV) 
  c1h = 0.000595747330615214     (/mV2) 
  d1h = 6.679127175885362e-06     (/mV3) 
  b2h = 0.10297542273524068     (/mV) 
  c2h = -0.001170310487510578     (/mV2) 
  d2h = 4.178083346880952e-06     (/mV3) 

  am = 0.09667136880137901     (/mV) 
  bm = -1.7473980140297374     (1) 
  vhm = -26.21286848058977     (mV) 
  Am = 0.23594832234835186     (/ms) 
  b1m = 0.04511591652368361     (/mV) 
  c1m = 0.0005527490983658916     (/mV2) 
  d1m = 3.0439610442900687e-06     (/mV3) 
  b2m = 0.027265009638420417     (/mV) 
  c2m = -0.00011870102107074166     (/mV2) 
  d2m = 1.291371606362106e-07     (/mV3) 
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