NEURON
{
  SUFFIX hh3 
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

  as = -0.33330455432551864     (/mV) 
  bs = 14.66552839388718     (1) 
  vhs = -29.983557998936554     (mV) 
  As = 50.65081092236301     (/ms) 
  b1s = -0.8580728174912685     (/mV) 
  c1s = 0.014265173511123094     (/mV2) 
  d1s = -6.198519195724137e-05     (/mV3) 
  b2s = -0.00025517399254494655     (/mV) 
  c2s = -6.567226857447942e-06     (/mV2) 
  d2s = -5.053103270966513e-08     (/mV3) 

  an = 0.33333266854126953     (/mV) 
  bn = -13.333320004057564     (1) 
  vhn = -27.421251637288787     (mV) 
  An = 2.0050991443535255     (/ms) 
  b1n = -0.01982222326666949     (/mV) 
  c1n = 0.000189538707568167     (/mV2) 
  d1n = -5.98839502434762e-07     (/mV3) 
  b2n = -0.020219567122600552     (/mV) 
  c2n = -0.00022492437228936933     (/mV2) 
  d2n = -9.68524417102729e-07     (/mV3) 

  ah = -0.33333590400249996     (/mV) 
  bh = 15.000118338934342     (1) 
  vhh = -100.0066655307623     (mV) 
  Ah = 1.0261244199713506     (/ms) 
  b1h = -0.0006178986641288609     (/mV) 
  c1h = 7.122848422957631e-06     (/mV2) 
  d1h = -2.1262349476931024e-08     (/mV3) 
  b2h = 0.0006202463204926525     (/mV) 
  c2h = -7.118845521992347e-06     (/mV2) 
  d2h = 2.124893784469861e-08     (/mV3) 

  am = 0.3333326685411616     (/mV) 
  bm = -13.333320004053096     (1) 
  vhm = -99.99990617158738     (mV) 
  Am = 0.12417328555583174     (/ms) 
  b1m = -0.0032181918603906203     (/mV) 
  c1m = 2.848836763380512e-05     (/mV2) 
  d1m = -9.378137993750374e-08     (/mV3) 
  b2m = 0.003218368915513632     (/mV) 
  c2m = -2.848853187436902e-05     (/mV2) 
  d2m = 9.378162491566885e-08     (/mV3) 

  an2 = 0.33333236619760814     (/mV) 
  bn2 = -13.333307708523094     (1) 
  vhn2 = -101.73687784226281     (mV) 
  An2 = 20.02484414677024     (/ms) 
  b1n2 = -1.7143310587088646e-05     (/mV) 
  c1n2 = 2.8639719244920566e-07     (/mV2) 
  d1n2 = -8.898346218529426e-10     (/mV3) 
  b2n2 = 3.1165421058779854e-05     (/mV) 
  c2n2 = -2.740722274010578e-07     (/mV2) 
  d2n2 = 7.954353039920924e-10     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  sInf 
  sTau 
  nInf 
  nTau 
  hInf 
  hTau 
  mInf 
  mTau 
  n2Inf 
  n2Tau 
}

STATE
{
  s
  n
  h
  m
  n2
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*s*n*n*h*m*n2
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  s' = (sInf - s) / sTau 
  n' = (nInf - n) / nTau 
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
  n2' = (n2Inf - n2) / n2Tau 
}

INITIAL
{
  rates(v)
  s = sInf 
  n = nInf 
  h = hInf 
  m = mInf 
  n2 = n2Inf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    sInf = 1/(1 + exp(-as*v + bs)) 
    sTau = As / ( exp(-(b1s*(v-vhs) + c1s*(v-vhs)^2 + d1s*(v-vhs)^3)) + exp((b2s*(v-vhs) + c2s*(v-vhs)^2 + d2s*(v-vhs)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    n2Inf = 1/(1 + exp(-an2*v + bn2)) 
    n2Tau = An2 / ( exp(-(b1n2*(v-vhn2) + c1n2*(v-vhn2)^2 + d1n2*(v-vhn2)^3)) + exp((b2n2*(v-vhn2) + c2n2*(v-vhn2)^2 + d2n2*(v-vhn2)^3)) ) 


  UNITSON
}