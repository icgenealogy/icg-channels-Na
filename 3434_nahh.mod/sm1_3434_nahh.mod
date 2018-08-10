NEURON
{
  SUFFIX nahh 
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

  ah = -0.14147105002096064     (/mV) 
  bh = 8.793782930948273     (1) 
  vhh = -67.24919895329549     (mV) 
  Ah = 0.605848962392589     (/ms) 
  b1h = 0.0978745815171717     (/mV) 
  c1h = 0.0024844085699179267     (/mV2) 
  d1h = 4.02859779634958e-05     (/mV3) 
  b2h = 0.0902248928298633     (/mV) 
  c2h = -0.0009023835689362094     (/mV2) 
  d2h = 2.7955955298804075e-06     (/mV3) 

  am = 0.10469485502915314     (/mV) 
  bm = -4.143238328779592     (1) 
  vhm = -51.65132491119346     (mV) 
  Am = 0.30602762481272044     (/ms) 
  b1m = -2.1089395569822913     (/mV) 
  c1m = 0.054477599325860326     (/mV2) 
  d1m = -0.0002738360062508367     (/mV3) 
  b2m = 0.007134402216433617     (/mV) 
  c2m = 0.0005595714814640929     (/mV2) 
  d2m = -3.7157451973462825e-06     (/mV3) 
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