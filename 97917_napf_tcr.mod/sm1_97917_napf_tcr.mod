NEURON
{
  SUFFIX napf_tcr 
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

  am = 0.10020930864753931     (/mV) 
  bm = -4.516508702612368     (1) 
  vhm = -40.62758112600384     (mV) 
  Am = 0.2644533337335219     (/ms) 
  b1m = 0.1761557992163961     (/mV) 
  c1m = 0.0047919863621385525     (/mV2) 
  d1m = 4.2015050793673104e-05     (/mV3) 
  b2m = 0.09873416654149732     (/mV) 
  c2m = -0.0012715496839636052     (/mV2) 
  d2m = 4.933968324938107e-06     (/mV3) 
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