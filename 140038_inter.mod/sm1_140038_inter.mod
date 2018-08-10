NEURON
{
  SUFFIX inter 
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

  ah = -0.12497279760970036     (/mV) 
  bh = 9.060725920870672     (1) 
  vhh = -105.13527215726592     (mV) 
  Ah = 308.1011474508143     (/ms) 
  b1h = 3.4373501121204293     (/mV) 
  c1h = 0.005052492006867123     (/mV2) 
  d1h = -1.319140654328158e-05     (/mV3) 
  b2h = 0.06881261731262228     (/mV) 
  c2h = 9.104529659751428e-07     (/mV2) 
  d2h = -1.6786459413196065e-08     (/mV3) 

  am = 0.10469571004718335     (/mV) 
  bm = -2.572840130698536     (1) 
  vhm = -106.73518609055445     (mV) 
  Am = 635.6035360246383     (/ms) 
  b1m = -0.11433798146839014     (/mV) 
  c1m = 0.0010526915717136945     (/mV2) 
  d1m = -7.209668700977682e-06     (/mV3) 
  b2m = -0.03464912759104434     (/mV) 
  c2m = 0.003349513723378923     (/mV2) 
  d2m = -0.00017548120372860745     (/mV3) 
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