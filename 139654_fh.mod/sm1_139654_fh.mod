NEURON
{
  SUFFIX fh 
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

  ah = -0.21488594521610302     (/mV) 
  bh = 13.481398530161693     (1) 
  vhh = -60.770584056182784     (mV) 
  Ah = 1.5374891297927484     (/ms) 
  b1h = -0.12369579620554436     (/mV) 
  c1h = 0.0012726012860188101     (/mV2) 
  d1h = -4.072733277901064e-06     (/mV3) 
  b2h = -0.10638945790944405     (/mV) 
  c2h = -0.0006336304039842476     (/mV2) 
  d2h = 3.979873533435403e-06     (/mV3) 

  am = 0.10378618890475713     (/mV) 
  bm = -3.4864173626099793     (1) 
  vhm = -22.727540022391054     (mV) 
  Am = 0.15676181132426653     (/ms) 
  b1m = 0.9934994276495539     (/mV) 
  c1m = 0.035003674776432896     (/mV2) 
  d1m = 0.0003184326914772344     (/mV3) 
  b2m = 0.08191633940439397     (/mV) 
  c2m = -0.001255426582358719     (/mV2) 
  d2m = 5.625447029621141e-06     (/mV3) 
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
  g = gbar*h*m*m
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