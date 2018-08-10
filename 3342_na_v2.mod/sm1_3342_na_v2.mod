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

  ah = -0.3778585132406626     (/mV) 
  bh = 18.12074427198179     (1) 
  vhh = -75.13299851985154     (mV) 
  Ah = 13.908004817788235     (/ms) 
  b1h = -0.038441398809107394     (/mV) 
  c1h = -4.372197825589875e-05     (/mV2) 
  d1h = 2.2127728697709026e-07     (/mV3) 
  b2h = -0.05474988281837369     (/mV) 
  c2h = -8.370096904057142e-05     (/mV2) 
  d2h = -6.558391335908499e-07     (/mV3) 

  am = 0.14799640098842431     (/mV) 
  bm = -4.622722396281253     (1) 
  vhm = 19.73134959661117     (mV) 
  Am = 0.23929383266697268     (/ms) 
  b1m = 0.018647643941153636     (/mV) 
  c1m = 2.522728781530843e-05     (/mV2) 
  d1m = 1.1050177339510118e-07     (/mV3) 
  b2m = 3.300748181664463     (/mV) 
  c2m = -0.02144125022772405     (/mV2) 
  d2m = 9.221147941194394e-05     (/mV3) 
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