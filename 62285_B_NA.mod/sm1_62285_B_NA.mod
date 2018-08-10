NEURON
{
  SUFFIX B_Na 
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

  ah = -0.11111099124128208     (/mV) 
  bh = 8.33332565666116     (1) 
  vhh = -54.60945790930994     (mV) 
  Ah = 1.1103889497304753     (/ms) 
  b1h = 0.01568552234138417     (/mV) 
  c1h = -2.9848471166521397e-05     (/mV2) 
  d1h = 5.6017898590865554e-08     (/mV3) 
  b2h = 0.23932659175367074     (/mV) 
  c2h = -0.004850382291066722     (/mV2) 
  d2h = 3.318838197431547e-05     (/mV3) 

  am = 0.11111143809322897     (/mV) 
  bm = -4.050418853574698     (1) 
  vhm = -51.70257476148035     (mV) 
  Am = 0.8381767138140278     (/ms) 
  b1m = -0.017796241597300054     (/mV) 
  c1m = -0.0010957202056745298     (/mV2) 
  d1m = -1.0813098729465527e-05     (/mV3) 
  b2m = 0.27574678278339815     (/mV) 
  c2m = -0.008675617963373846     (/mV2) 
  d2m = 9.28649637921784e-05     (/mV3) 
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