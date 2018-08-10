NEURON
{
  SUFFIX ls 
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

  ah = -0.14633513003672255     (/mV) 
  bh = 8.20199576483212     (1) 
  vhh = -104.6865186233761     (mV) 
  Ah = 3154.9629755407786     (/ms) 
  b1h = -0.06437839993983986     (/mV) 
  c1h = 0.00035959793933350363     (/mV2) 
  d1h = -6.76604725694687e-07     (/mV3) 
  b2h = 0.009843441596018541     (/mV) 
  c2h = -0.000321632962020551     (/mV2) 
  d2h = -5.871209541730925e-06     (/mV3) 

  am = 0.2173674309064108     (/mV) 
  bm = -11.264727618083562     (1) 
  vhm = -133.6052161653054     (mV) 
  Am = 4237.134759574144     (/ms) 
  b1m = 0.43409514766268226     (/mV) 
  c1m = -3.5829511615188784e-05     (/mV2) 
  d1m = -2.0667245860180058e-05     (/mV3) 
  b2m = 0.0323392136759589     (/mV) 
  c2m = 9.836047507085057e-05     (/mV2) 
  d2m = -8.375652706808423e-07     (/mV3) 
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
  g = gbar*h*m
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