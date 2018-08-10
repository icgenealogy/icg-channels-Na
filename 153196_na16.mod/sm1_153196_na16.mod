NEURON
{
  SUFFIX na16 
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

  ah = -0.16110831096426437     (/mV) 
  bh = 10.79424793177341     (1) 
  vhh = -61.80136064804162     (mV) 
  Ah = 19.212436790059513     (/ms) 
  b1h = 0.14671921383930492     (/mV) 
  c1h = 0.0029002287811235937     (/mV2) 
  d1h = 2.6309971770612447e-05     (/mV3) 
  b2h = 0.14185716626520772     (/mV) 
  c2h = -0.0014991413898000196     (/mV2) 
  d2h = 5.4164800880267446e-06     (/mV3) 

  am = 0.1664714364562524     (/mV) 
  bm = -6.709235504777538     (1) 
  vhm = -46.064781439286016     (mV) 
  Am = 0.3435596479514336     (/ms) 
  b1m = -0.06027300346545766     (/mV) 
  c1m = 0.0004582254358012716     (/mV2) 
  d1m = -1.24993643629886e-06     (/mV3) 
  b2m = -0.09084524781128508     (/mV) 
  c2m = -0.0015690535596363105     (/mV2) 
  d2m = -1.0444261572908726e-05     (/mV3) 
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