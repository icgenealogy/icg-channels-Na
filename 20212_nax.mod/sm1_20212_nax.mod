NEURON
{
  SUFFIX nax 
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

  ah = -0.24999955331037577     (/mV) 
  bh = 12.499975673542856     (1) 
  vhh = -48.588211938323816     (mV) 
  Ah = 15.020379746520636     (/ms) 
  b1h = -0.19841291761395918     (/mV) 
  c1h = 0.002760874137099447     (/mV2) 
  d1h = -1.103018246587544e-05     (/mV3) 
  b2h = -0.22579872724729372     (/mV) 
  c2h = -0.0071527961627621435     (/mV2) 
  d2h = -7.864845094713633e-05     (/mV3) 

  am = 0.13888872579100658     (/mV) 
  bm = -5.337837931299445     (1) 
  vhm = -33.22076321875721     (mV) 
  Am = 0.24453680855779153     (/ms) 
  b1m = 0.02872048633845459     (/mV) 
  c1m = -4.8858328561324556e-05     (/mV2) 
  d1m = -2.031721975744898e-06     (/mV3) 
  b2m = 0.07908667846228963     (/mV) 
  c2m = -0.0007873040389808289     (/mV2) 
  d2m = 2.2065926644326795e-06     (/mV3) 
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