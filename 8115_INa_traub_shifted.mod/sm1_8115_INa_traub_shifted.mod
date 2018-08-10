NEURON
{
  SUFFIX inaT 
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

  ah = -0.24748563538836624     (/mV) 
  bh = 10.7209362968196     (1) 
  vhh = -44.33441721768196     (mV) 
  Ah = 10.59130731983525     (/ms) 
  b1h = -0.1952615384335411     (/mV) 
  c1h = 0.002645137300215533     (/mV2) 
  d1h = -1.0523408874596936e-05     (/mV3) 
  b2h = -0.07186501058407652     (/mV) 
  c2h = -0.0006943224075344804     (/mV2) 
  d2h = -8.242585944370642e-06     (/mV3) 

  am = 0.1349151345677765     (/mV) 
  bm = -5.265108161490034     (1) 
  vhm = -56.099715947252264     (mV) 
  Am = 0.2150796270727938     (/ms) 
  b1m = -0.024626431575801987     (/mV) 
  c1m = 5.91945035503155e-05     (/mV2) 
  d1m = 1.1204646563818896e-07     (/mV3) 
  b2m = -0.05748364465052122     (/mV) 
  c2m = -0.0009276183034911499     (/mV2) 
  d2m = -7.124186164279236e-06     (/mV3) 
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