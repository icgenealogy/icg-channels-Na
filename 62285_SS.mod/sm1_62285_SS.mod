NEURON
{
  SUFFIX SS 
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

  ah = -0.11112906274742157     (/mV) 
  bh = 8.334754661760707     (1) 
  vhh = -67.24012012167034     (mV) 
  Ah = 2.1202154229151016     (/ms) 
  b1h = 0.15488042229746396     (/mV) 
  c1h = 0.003568304841893754     (/mV2) 
  d1h = 4.0384291565769404e-05     (/mV3) 
  b2h = 0.13685022674071007     (/mV) 
  c2h = -0.001387304434176121     (/mV2) 
  d2h = 4.627300408581634e-06     (/mV3) 

  am = 0.11136543029324769     (/mV) 
  bm = -5.399065611718965     (1) 
  vhm = -93.2955139609444     (mV) 
  Am = 0.022031573820976716     (/ms) 
  b1m = 6.251748028301209     (/mV) 
  c1m = -0.07288228664307615     (/mV2) 
  d1m = -0.0002006200900277203     (/mV3) 
  b2m = 0.08750603473591399     (/mV) 
  c2m = -0.0036562545475132863     (/mV2) 
  d2m = 3.814730502612164e-05     (/mV3) 
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
  g = gbar*m*m*m
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