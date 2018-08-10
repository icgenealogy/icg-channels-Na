NEURON
{
  SUFFIX ichan 
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

  ah = -0.14591410005661412     (/mV) 
  bh = 6.501011162060365     (1) 
  vhh = -41.40323570300264     (mV) 
  Ah = 0.33019324495956115     (/ms) 
  b1h = -0.1162508024064084     (/mV) 
  c1h = 0.0014201118016707252     (/mV2) 
  d1h = -5.308050703987045e-06     (/mV3) 
  b2h = -0.03762720300873873     (/mV) 
  c2h = 0.00013196584860337114     (/mV2) 
  d2h = 5.746364882790402e-07     (/mV3) 

  anf = 0.12346337026155864     (/mV) 
  bnf = -3.272114190296453     (1) 
  vhnf = -25.43918915703582     (mV) 
  Anf = 0.18398522406760032     (/ms) 
  b1nf = -0.10653614066809221     (/mV) 
  c1nf = 0.001976160275610243     (/mV2) 
  d1nf = -1.5247984159515117e-05     (/mV3) 
  b2nf = -0.017125092194956566     (/mV) 
  c2nf = 0.0001663843245489111     (/mV2) 
  d2nf = 1.5887675537941744e-06     (/mV3) 

  am = 0.12492610805973998     (/mV) 
  bm = -3.2481657713500347     (1) 
  vhm = -52.047471580720085     (mV) 
  Am = 0.3735680157786723     (/ms) 
  b1m = 0.12417301725299253     (/mV) 
  c1m = 0.007391805124455478     (/mV2) 
  d1m = 0.0002372237366635402     (/mV3) 
  b2m = 0.07050906304307188     (/mV) 
  c2m = -0.00023600507911482438     (/mV2) 
  d2m = 2.2684809943570852e-06     (/mV3) 
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
}

STATE
{
  h
  nf
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*nf*nf*nf*nf*m*m*m
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  nf' = (nfInf - nf) / nfTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  nf = nfInf 
  m = mInf 
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


  UNITSON
}