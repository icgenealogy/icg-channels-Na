NEURON
{
  SUFFIX ichan2 
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

  ah = -0.1459140719334338     (/mV) 
  bh = 6.938754446399779     (1) 
  vhh = -44.14439882399019     (mV) 
  Ah = 9.274102007665604     (/ms) 
  b1h = 0.0350480540206678     (/mV) 
  c1h = -0.00020186654826211947     (/mV2) 
  d1h = -2.9814828653851547e-07     (/mV3) 
  b2h = 0.12133901311759075     (/mV) 
  c2h = -0.0013401734193643656     (/mV2) 
  d2h = 4.642346708084336e-06     (/mV3) 

  anf = 0.12346211175367248     (/mV) 
  bnf = -3.2720944914787213     (1) 
  vhnf = -29.287532175333446     (mV) 
  Anf = 6.049988061920639     (/ms) 
  b1nf = -0.09923136398212175     (/mV) 
  c1nf = 0.0010447132425556581     (/mV2) 
  d1nf = -4.033429265485314e-06     (/mV3) 
  b2nf = -0.034544446630696994     (/mV) 
  c2nf = -0.0002266041077945346     (/mV2) 
  d2nf = -1.6586135438922715e-06     (/mV3) 

  am = 0.1249260871399721     (/mV) 
  bm = -3.6229475480161115     (1) 
  vhm = -34.86329658157132     (mV) 
  Am = 0.2292939553968912     (/ms) 
  b1m = -0.03773332046669612     (/mV) 
  c1m = 0.00027307118759936817     (/mV2) 
  d1m = -8.110744669882462e-07     (/mV3) 
  b2m = -0.05484467321830793     (/mV) 
  c2m = -0.0008159509976556636     (/mV2) 
  d2m = -5.772887994467185e-06     (/mV3) 

  ans = 0.123464659297432     (/mV) 
  bns = -4.753692515053472     (1) 
  vhns = -44.61166618361986     (mV) 
  Ans = 16.521973921582443     (/ms) 
  b1ns = -0.09333960222526408     (/mV) 
  c1ns = 0.0008726941714666216     (/mV2) 
  d1ns = -2.9942816702364235e-06     (/mV3) 
  b2ns = -0.05014906001742162     (/mV) 
  c2ns = -0.0006953571045396485     (/mV2) 
  d2ns = -6.174148549550419e-06     (/mV3) 
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
  nsInf 
  nsTau 
}

STATE
{
  h
  nf
  m
  ns
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*nf*nf*nf*nf*m*m*m*ns*ns*ns*ns
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  nf' = (nfInf - nf) / nfTau 
  m' = (mInf - m) / mTau 
  ns' = (nsInf - ns) / nsTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  nf = nfInf 
  m = mInf 
  ns = nsInf 
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

    nsInf = 1/(1 + exp(-ans*v + bns)) 
    nsTau = Ans / ( exp(-(b1ns*(v-vhns) + c1ns*(v-vhns)^2 + d1ns*(v-vhns)^3)) + exp((b2ns*(v-vhns) + c2ns*(v-vhns)^2 + d2ns*(v-vhns)^3)) ) 


  UNITSON
}