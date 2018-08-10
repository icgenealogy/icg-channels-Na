NEURON
{
  SUFFIX Nap_Et2 
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

  ah = -0.0995434041330054     (/mV) 
  bh = 4.863181065066581     (1) 
  vhh = -7.648620487165871     (mV) 
  Ah = 1638.7638809775856     (/ms) 
  b1h = -0.01460647864843899     (/mV) 
  c1h = 0.0003600571828761808     (/mV2) 
  d1h = -2.5740115267846447e-06     (/mV3) 
  b2h = 0.014594579382435836     (/mV) 
  c2h = 9.087383286135682e-05     (/mV2) 
  d2h = -1.0121424647122892e-06     (/mV3) 

  am = 0.2173909738721496     (/mV) 
  bm = -11.43476215929012     (1) 
  vhm = -44.606765803021126     (mV) 
  Am = 2.2403003330113203     (/ms) 
  b1m = -0.06721779053752712     (/mV) 
  c1m = 0.0005712644404824684     (/mV2) 
  d1m = -1.8446385935695236e-06     (/mV3) 
  b2m = -0.09146690270851784     (/mV) 
  c2m = -0.0016635376292980521     (/mV2) 
  d2m = -1.2878234398668241e-05     (/mV3) 
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