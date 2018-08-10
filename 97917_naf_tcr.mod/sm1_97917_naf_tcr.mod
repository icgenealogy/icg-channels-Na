NEURON
{
  SUFFIX naf_tcr 
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

  ah = -0.0934613976083913     (/mV) 
  bh = 5.224552506272413     (1) 
  vhh = -18.245747381490528     (mV) 
  Ah = 0.8195644715195853     (/ms) 
  b1h = -0.05375033346721801     (/mV) 
  c1h = 0.0007313540302715778     (/mV2) 
  d1h = -2.756695226245081e-06     (/mV3) 
  b2h = 0.009488797522678285     (/mV) 
  c2h = 3.0319608362962117e-05     (/mV2) 
  d2h = -2.2497661039709964e-07     (/mV3) 

  am = 0.10000303276828658     (/mV) 
  bm = -3.800327590614767     (1) 
  vhm = -32.69829768967778     (mV) 
  Am = 0.27918667335852015     (/ms) 
  b1m = -0.10801516303329504     (/mV) 
  c1m = 0.001482393474974907     (/mV2) 
  d1m = -6.122628781223566e-06     (/mV3) 
  b2m = -0.1720688010814768     (/mV) 
  c2m = -0.004269430211062446     (/mV2) 
  d2m = -3.3460700626402436e-05     (/mV3) 
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