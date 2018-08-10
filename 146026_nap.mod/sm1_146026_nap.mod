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

  am = 0.09999702286293913     (/mV) 
  bm = -4.799854502105199     (1) 
  vhm = -43.511172991933016     (mV) 
  Am = 0.27243928067046747     (/ms) 
  b1m = -0.10001659403322335     (/mV) 
  c1m = 0.0012773770870492926     (/mV2) 
  d1m = -4.898864871791468e-06     (/mV3) 
  b2m = -0.18955130352077804     (/mV) 
  c2m = -0.005506235709333868     (/mV2) 
  d2m = -5.1534567822905996e-05     (/mV3) 
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