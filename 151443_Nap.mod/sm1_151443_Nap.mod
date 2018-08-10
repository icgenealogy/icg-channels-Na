NEURON
{
  SUFFIX Nap 
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

  ap = 0.21151310610601712     (/mV) 
  bp = -8.752444244315312     (1) 
  vhp = -41.37157583050737     (mV) 
  Ap = 26.867532959825766     (/ms) 
  b1p = -0.1374397537618778     (/mV) 
  c1p = 0.0016021929502571409     (/mV2) 
  d1p = -6.478012798471502e-06     (/mV3) 
  b2p = -0.06893020982340717     (/mV) 
  c2p = -0.001071819786805057     (/mV2) 
  d2p = -7.4938151301885036e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  pInf 
  pTau 
}

STATE
{
  p
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*p*p*p
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  p' = (pInf - p) / pTau 
}

INITIAL
{
  rates(v)
  p = pInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    pInf = 1/(1 + exp(-ap*v + bp)) 
    pTau = Ap / ( exp(-(b1p*(v-vhp) + c1p*(v-vhp)^2 + d1p*(v-vhp)^3)) + exp((b2p*(v-vhp) + c2p*(v-vhp)^2 + d2p*(v-vhp)^3)) ) 


  UNITSON
}