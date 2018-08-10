NEURON
{
  SUFFIX navn 
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

  ah = -0.49999876513063796     (/mV) 
  bh = 29.999921580054277     (1) 
  vhh = -85.98108029949947     (mV) 
  Ah = 9.232585465994621     (/ms) 
  b1h = -0.028314299087998807     (/mV) 
  c1h = -0.00030947926917528285     (/mV2) 
  d1h = 1.88562824787771e-06     (/mV3) 
  b2h = -0.051563422453756634     (/mV) 
  c2h = -5.354162432701965e-07     (/mV2) 
  d2h = 2.1619160364831994e-06     (/mV3) 

  am = 0.12499989351294442     (/mV) 
  bm = -5.249990242296119     (1) 
  vhm = -83.33533233138617     (mV) 
  Am = 1.1177529264569073     (/ms) 
  b1m = 0.004755137452532298     (/mV) 
  c1m = 0.00020926086114259846     (/mV2) 
  d1m = -2.404316880283243e-06     (/mV3) 
  b2m = -0.00476098219307903     (/mV) 
  c2m = 0.00042911476098222187     (/mV2) 
  d2m = -1.5006021503170064e-06     (/mV3) 
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