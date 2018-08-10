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

  am = 0.3333326685412292     (/mV) 
  bm = -13.333320004058379     (1) 
  vhm = -99.99990617158738     (mV) 
  Am = 0.12417328555583174     (/ms) 
  b1m = -0.0032181918603906203     (/mV) 
  c1m = 2.848836763380512e-05     (/mV2) 
  d1m = -9.378137993750374e-08     (/mV3) 
  b2m = 0.003218368915513632     (/mV) 
  c2m = -2.848853187436902e-05     (/mV2) 
  d2m = 9.378162491566885e-08     (/mV3) 

  ah = -0.33333590400244206     (/mV) 
  bh = 15.000118338931943     (1) 
  vhh = -100.0066655307623     (mV) 
  Ah = 1.0261244199713506     (/ms) 
  b1h = -0.0006178986641288609     (/mV) 
  c1h = 7.122848422957631e-06     (/mV2) 
  d1h = -2.1262349476931024e-08     (/mV3) 
  b2h = 0.0006202463204926525     (/mV) 
  c2h = -7.118845521992347e-06     (/mV2) 
  d2h = 2.124893784469861e-08     (/mV3) 

  an = 0.3333326685412009     (/mV) 
  bn = -13.333320004057258     (1) 
  vhn = -27.421251637288787     (mV) 
  An = 2.0050991443535255     (/ms) 
  b1n = -0.01982222326666949     (/mV) 
  c1n = 0.000189538707568167     (/mV2) 
  d1n = -5.98839502434762e-07     (/mV3) 
  b2n = -0.020219567122600552     (/mV) 
  c2n = -0.00022492437228936933     (/mV2) 
  d2n = -9.68524417102729e-07     (/mV3) 

  as = -0.3333045543225761     (/mV) 
  bs = 14.66552839375724     (1) 
  vhs = -29.983557998936554     (mV) 
  As = 50.65081092236301     (/ms) 
  b1s = -0.8580728174912685     (/mV) 
  c1s = 0.014265173511123094     (/mV2) 
  d1s = -6.198519195724137e-05     (/mV3) 
  b2s = -0.00025517399254494655     (/mV) 
  c2s = -6.567226857447942e-06     (/mV2) 
  d2s = -5.053103270966513e-08     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
  hInf 
  hTau 
  nInf 
  nTau 
  sInf 
  sTau 
}

STATE
{
  m
  h
  n
  s
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*h*n*s
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
  s' = (sInf - s) / sTau 
}

INITIAL
{
  rates(v)
  m = mInf 
  h = hInf 
  n = nInf 
  s = sInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    sInf = 1/(1 + exp(-as*v + bs)) 
    sTau = As / ( exp(-(b1s*(v-vhs) + c1s*(v-vhs)^2 + d1s*(v-vhs)^3)) + exp((b2s*(v-vhs) + c2s*(v-vhs)^2 + d2s*(v-vhs)^3)) ) 


  UNITSON
}