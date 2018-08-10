NEURON
{
  SUFFIX NaPm 
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

  am = 0.3225797577278361     (/mV) 
  bm = -15.4193054807441     (1) 
  vhm = -61.286462867376045     (mV) 
  Am = 0.5200000003979275     (/ms) 
  b1m = -2.7276648320281297e-09     (/mV) 
  c1m = 5.546166329054749e-10     (/mV2) 
  d1m = 4.5278765934308743e-11     (/mV3) 
  b2m = -2.731955978592326e-09     (/mV) 
  c2m = 5.566312776564243e-10     (/mV2) 
  d2m = 4.5257903285997825e-11     (/mV3) 
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