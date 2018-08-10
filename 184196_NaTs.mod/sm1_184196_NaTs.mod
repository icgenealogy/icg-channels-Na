NEURON
{
  SUFFIX NaTs 
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

  ah = -0.16666639513868173     (/mV) 
  bh = 10.999981729139735     (1) 
  vhh = -72.14834471373017     (mV) 
  Ah = 3.2030942121036183     (/ms) 
  b1h = -0.05465651989693992     (/mV) 
  c1h = 0.0003852951648562962     (/mV2) 
  d1h = -1.0421627026165824e-06     (/mV3) 
  b2h = -0.1063958329211918     (/mV) 
  c2h = -0.002693544140869408     (/mV2) 
  d2h = -3.527872405214354e-05     (/mV3) 

  am = 0.1666663046640272     (/mV) 
  bm = -7.05037082727696     (1) 
  vhm = -47.87584907206504     (mV) 
  Am = 0.3407367375206721     (/ms) 
  b1m = -0.061775450223582606     (/mV) 
  c1m = 0.0004954373349331031     (/mV2) 
  d1m = -1.4367967957071703e-06     (/mV3) 
  b2m = -0.09424795447906739     (/mV) 
  c2m = -0.0019582675609606573     (/mV2) 
  d2m = -1.6845422454067004e-05     (/mV3) 
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