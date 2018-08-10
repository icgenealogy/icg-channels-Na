NEURON
{
  SUFFIX nahh 
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

  ah = -0.19999335947537358     (/mV) 
  bh = 10.399692699101685     (1) 
  vhh = -52.32373073493683     (mV) 
  Ah = 46.560000024334066     (/ms) 
  b1h = -8.98305121584614e-09     (/mV) 
  c1h = 4.164816433975987e-10     (/mV2) 
  d1h = -5.116089362353202e-11     (/mV3) 
  b2h = -8.957180608921424e-09     (/mV) 
  c2h = 4.176064606616427e-10     (/mV2) 
  d2h = -5.117652255335396e-11     (/mV3) 

  am = 0.142856902000728     (/mV) 
  bm = -3.5714208560325567     (1) 
  vhm = -36.99196464425484     (mV) 
  Am = 2.800000000743788     (/ms) 
  b1m = 4.625500287451176e-08     (/mV) 
  c1m = -4.023922000925017e-09     (/mV2) 
  d1m = -5.277806292483675e-11     (/mV3) 
  b2m = 4.632629803598757e-08     (/mV) 
  c2m = -4.0233754666223495e-09     (/mV2) 
  d2m = -5.280015405800723e-11     (/mV3) 
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