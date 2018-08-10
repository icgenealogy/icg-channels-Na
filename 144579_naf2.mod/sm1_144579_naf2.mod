NEURON
{
  SUFFIX naf2 
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

  ah = -0.14924589052379314     (/mV) 
  bh = 8.701043147920268     (1) 
  vhh = -48.39997339140547     (mV) 
  Ah = 1.9996304072026183     (/ms) 
  b1h = -0.061324735483537736     (/mV) 
  c1h = 0.0005730164987941781     (/mV2) 
  d1h = -1.7346981732337213e-06     (/mV3) 
  b2h = -0.028101633772168805     (/mV) 
  c2h = -0.0006685191893632568     (/mV2) 
  d2h = -5.307838580870396e-06     (/mV3) 

  am = 0.09999623422367836     (/mV) 
  bm = -3.4498974363660744     (1) 
  vhm = -29.19148795294483     (mV) 
  Am = 0.2605831547718765     (/ms) 
  b1m = -0.10413713940978582     (/mV) 
  c1m = 0.00144313985570351     (/mV2) 
  d1m = -6.054551274191392e-06     (/mV3) 
  b2m = -0.17876446086458514     (/mV) 
  c2m = -0.004048820287201607     (/mV2) 
  d2m = -2.9304528448560162e-05     (/mV3) 
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