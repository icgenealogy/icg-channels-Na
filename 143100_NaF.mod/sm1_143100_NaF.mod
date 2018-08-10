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

  ah = -0.3571418510037951     (/mV) 
  bh = 17.142799519485724     (1) 
  vhh = -47.38819861308023     (mV) 
  Ah = 4.399424020111526     (/ms) 
  b1h = -0.13571630154857098     (/mV) 
  c1h = 0.0017511703285467025     (/mV2) 
  d1h = -6.689001033348548e-06     (/mV3) 
  b2h = -0.1709335248027823     (/mV) 
  c2h = -0.0041286868361305774     (/mV2) 
  d2h = -3.688590300732678e-05     (/mV3) 

  am = 0.1999996351905595     (/mV) 
  bm = -7.799977037152167     (1) 
  vhm = -59.14115325919962     (mV) 
  Am = 0.060000000117850345     (/ms) 
  b1m = 6.973269594815698e-12     (/mV) 
  c1m = -1.065700646529233e-08     (/mV2) 
  d1m = -3.612990854852529e-11     (/mV3) 
  b2m = 3.1251979726860316e-11     (/mV) 
  c2m = -1.0651333960969128e-08     (/mV2) 
  d2m = -3.620610875074159e-11     (/mV3) 

  as = -0.10413146855837015     (/mV) 
  bs = 3.738328220739713     (1) 
  vhs = -37.024482913595975     (mV) 
  As = 924.9896468040095     (/ms) 
  b1s = -0.12080829727463238     (/mV) 
  c1s = 0.0008382909523276456     (/mV2) 
  d1s = -1.0287868225104064e-06     (/mV3) 
  b2s = -0.03987144088369168     (/mV) 
  c2s = 0.0002556743548545696     (/mV2) 
  d2s = 2.034626464291266e-06     (/mV3) 
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
  sInf 
  sTau 
}

STATE
{
  h
  m
  s
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m*m*m*s
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
  s' = (sInf - s) / sTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
  s = sInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    sInf = 1/(1 + exp(-as*v + bs)) 
    sTau = As / ( exp(-(b1s*(v-vhs) + c1s*(v-vhs)^2 + d1s*(v-vhs)^3)) + exp((b2s*(v-vhs) + c2s*(v-vhs)^2 + d2s*(v-vhs)^3)) ) 


  UNITSON
}