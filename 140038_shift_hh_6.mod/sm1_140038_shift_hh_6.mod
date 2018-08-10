NEURON
{
  SUFFIX Shh_6 
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

  am = 0.10469775952623359     (/mV) 
  bm = -3.5152231481261573     (1) 
  vhm = -29.505176786368846     (mV) 
  Am = 0.2746670109383073     (/ms) 
  b1m = -0.08424350613038739     (/mV) 
  c1m = 0.0001728450844386771     (/mV2) 
  d1m = 2.3387508068856885e-05     (/mV3) 
  b2m = 0.01874862005868038     (/mV) 
  c2m = 0.0003268723679908327     (/mV2) 
  d2m = -2.8774912058568333e-06     (/mV3) 

  ah = -0.14147695879399416     (/mV) 
  bh = 7.945293341304475     (1) 
  vhh = -60.697019443293385     (mV) 
  Ah = 0.6088694877064401     (/ms) 
  b1h = 0.09398566367465913     (/mV) 
  c1h = 0.0020833212948587923     (/mV2) 
  d1h = 2.8226954434335632e-05     (/mV3) 
  b2h = 0.09280256480646958     (/mV) 
  c2h = -0.0009566989140767863     (/mV2) 
  d2h = 3.0620634075900675e-06     (/mV3) 

  an = 0.05640514082345872     (/mV) 
  bn = -2.8790033405797906     (1) 
  vhn = -81.91485211966854     (mV) 
  An = 0.3991683192025798     (/ms) 
  b1n = 0.038359520638365445     (/mV) 
  c1n = 0.0007470300267753836     (/mV2) 
  d1n = 1.2940312722420579e-05     (/mV3) 
  b2n = 0.03343351304321064     (/mV) 
  c2n = -0.0001494471549557467     (/mV2) 
  d2n = 2.508970605067338e-07     (/mV3) 
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
}

STATE
{
  m
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*m*h*n*n*n*n
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  m = mInf 
  h = hInf 
  n = nInf 
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


  UNITSON
}