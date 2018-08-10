NEURON
{
  SUFFIX NaTa 
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

  ah = -0.1666663984724896     (/mV) 
  bh = 11.499981842919953     (1) 
  vhh = -75.18241245039847     (mV) 
  Ah = 3.1921767581636815     (/ms) 
  b1h = -0.05429823056130963     (/mV) 
  c1h = 0.0003790220121877964     (/mV2) 
  d1h = -1.011578415670937e-06     (/mV3) 
  b2h = -0.1071364016318212     (/mV) 
  c2h = -0.0028290330918182326     (/mV2) 
  d2h = -3.985961030163651e-05     (/mV3) 

  am = 0.1666663650636806     (/mV) 
  bm = -8.383705550795218     (1) 
  vhm = -57.20869362873134     (mV) 
  Am = 0.33664003906907053     (/ms) 
  b1m = -0.05668578416193969     (/mV) 
  c1m = 0.000393719972457748     (/mV2) 
  d1m = -9.49746549882083e-07     (/mV3) 
  b2m = -0.09840099080853992     (/mV) 
  c2m = -0.0021659000759686113     (/mV2) 
  d2m = -2.1140547508389585e-05     (/mV3) 
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