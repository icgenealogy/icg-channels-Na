NEURON
{
  SUFFIX NaF 
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

  ah = -0.16949106507729966     (/mV) 
  bh = 7.118618346220029     (1) 
  vhh = -53.931009061568034     (mV) 
  Ah = 3.129959846092265     (/ms) 
  b1h = 0.0737450227792797     (/mV) 
  c1h = 0.0016011776861207117     (/mV2) 
  d1h = 1.7179370928586294e-05     (/mV3) 
  b2h = 0.09792493875476938     (/mV) 
  c2h = -0.0010665467013391847     (/mV2) 
  d2h = 3.593628412225173e-06     (/mV3) 

  am = 0.1369859202713891     (/mV) 
  bm = -6.164362410723376     (1) 
  vhm = -32.278498645635494     (mV) 
  Am = 0.1557902320800909     (/ms) 
  b1m = 0.10163656876431282     (/mV) 
  c1m = 0.002045578145760441     (/mV2) 
  d1m = 1.3407412601780076e-05     (/mV3) 
  b2m = 0.05016868827984868     (/mV) 
  c2m = -0.0004967285163407388     (/mV2) 
  d2m = 1.6039662210736607e-06     (/mV3) 
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