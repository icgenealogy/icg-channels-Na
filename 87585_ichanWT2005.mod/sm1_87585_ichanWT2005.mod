NEURON
{
  SUFFIX ichanWT2005 
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

  ah = -0.14925398145774246     (/mV) 
  bh = 6.253739747434996     (1) 
  vhh = -20.38379637396055     (mV) 
  Ah = 0.687946975703534     (/ms) 
  b1h = -0.002360497604559656     (/mV) 
  c1h = 0.00012516892307277484     (/mV2) 
  d1h = 1.9016465163728615e-06     (/mV3) 
  b2h = 0.05759995708942248     (/mV) 
  c2h = -9.704164345734715e-05     (/mV2) 
  d2h = -4.4812552815326374e-07     (/mV3) 

  anf = 0.12346058389006537     (/mV) 
  bnf = -3.2719530863900435     (1) 
  vhnf = -27.896490343268063     (mV) 
  Anf = 0.20311760483298152     (/ms) 
  b1nf = -0.10185288138088808     (/mV) 
  c1nf = 0.001562821259464002     (/mV2) 
  d1nf = -1.0498271365786867e-05     (/mV3) 
  b2nf = -0.026903939694436538     (/mV) 
  c2nf = -2.3297863364141193e-05     (/mV2) 
  d2nf = 4.7326333748628455e-07     (/mV3) 

  am = 0.1850386242201919     (/mV) 
  bm = -5.070059041387868     (1) 
  vhm = -60.72644936075323     (mV) 
  Am = 1.3934707757086533     (/ms) 
  b1m = 1.444368712736585     (/mV) 
  c1m = -0.019909233224051616     (/mV2) 
  d1m = 3.768597494638218e-05     (/mV3) 
  b2m = 0.035155463403916905     (/mV) 
  c2m = -6.72698774740514e-05     (/mV2) 
  d2m = 9.963449244714852e-07     (/mV3) 

  as = -0.1498399449710049     (/mV) 
  bs = 6.904170073232951     (1) 
  vhs = -45.797543637096105     (mV) 
  As = 6850.526481848961     (/ms) 
  b1s = -0.021525732778871244     (/mV) 
  c1s = -0.00050946524723221     (/mV2) 
  d1s = 1.4913182731150106e-08     (/mV3) 
  b2s = 0.03189333819980185     (/mV) 
  c2s = 0.000510998819817051     (/mV2) 
  d2s = 6.847248194587999e-09     (/mV3) 
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
  nfInf 
  nfTau 
  mInf 
  mTau 
  sInf 
  sTau 
}

STATE
{
  h
  nf
  m
  s
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*nf*nf*nf*nf*m*m*m*s
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  nf' = (nfInf - nf) / nfTau 
  m' = (mInf - m) / mTau 
  s' = (sInf - s) / sTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  nf = nfInf 
  m = mInf 
  s = sInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nfInf = 1/(1 + exp(-anf*v + bnf)) 
    nfTau = Anf / ( exp(-(b1nf*(v-vhnf) + c1nf*(v-vhnf)^2 + d1nf*(v-vhnf)^3)) + exp((b2nf*(v-vhnf) + c2nf*(v-vhnf)^2 + d2nf*(v-vhnf)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    sInf = 1/(1 + exp(-as*v + bs)) 
    sTau = As / ( exp(-(b1s*(v-vhs) + c1s*(v-vhs)^2 + d1s*(v-vhs)^3)) + exp((b2s*(v-vhs) + c2s*(v-vhs)^2 + d2s*(v-vhs)^3)) ) 


  UNITSON
}