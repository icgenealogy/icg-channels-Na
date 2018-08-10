NEURON
{
  SUFFIX na3 
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

  ah = -0.16110658738082528     (/mV) 
  bh = 10.471933157031541     (1) 
  vhh = -66.95816941838565     (mV) 
  Ah = 19.211215841555223     (/ms) 
  b1h = 0.1502909566652062     (/mV) 
  c1h = 0.00317353638290402     (/mV2) 
  d1h = 3.2042051178312994e-05     (/mV3) 
  b2h = 0.1403002793404687     (/mV) 
  c2h = -0.0014547116975574357     (/mV2) 
  d2h = 5.141274352920912e-06     (/mV3) 

  am = 0.11105668402315941     (/mV) 
  bm = -4.2705568134419964     (1) 
  vhm = -45.94568255222068     (mV) 
  Am = 0.23406055292377392     (/ms) 
  b1m = -0.041725556219116465     (/mV) 
  c1m = 0.00023266386832819475     (/mV2) 
  d1m = -4.394037448025439e-07     (/mV3) 
  b2m = -0.06469408128758541     (/mV) 
  c2m = -0.0010169499671365643     (/mV2) 
  d2m = -7.518718858915409e-06     (/mV3) 
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