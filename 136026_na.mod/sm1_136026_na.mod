NEURON
{
  SUFFIX na 
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

  ah = -0.2670163139797174     (/mV) 
  bh = 12.865150952052097     (1) 
  vhh = -56.6459958837636     (mV) 
  Ah = 13.653320034086764     (/ms) 
  b1h = -0.1240181595779913     (/mV) 
  c1h = 0.00123277290906209     (/mV2) 
  d1h = -4.2809192702876395e-06     (/mV3) 
  b2h = -0.05721421457540868     (/mV) 
  c2h = -0.0010642271074129603     (/mV2) 
  d2h = -9.70650068867959e-06     (/mV3) 

  am = 0.10731206474496659     (/mV) 
  bm = -3.6124223595221454     (1) 
  vhm = -36.23932966405484     (mV) 
  Am = 0.10947376048738929     (/ms) 
  b1m = 0.064963092478985     (/mV) 
  c1m = 0.001233955545684285     (/mV2) 
  d1m = 1.0246048728622252e-05     (/mV3) 
  b2m = 0.017624245448545618     (/mV) 
  c2m = -8.151684736221975e-05     (/mV2) 
  d2m = 1.8082775422209968e-07     (/mV3) 
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