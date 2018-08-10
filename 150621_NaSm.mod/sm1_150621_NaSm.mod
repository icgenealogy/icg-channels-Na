NEURON
{
  SUFFIX NaSm 
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

  am = 0.10637917993997711     (/mV) 
  bm = -1.7021359141219576     (1) 
  vhm = -33.49823644536998     (mV) 
  Am = 147.2321444085874     (/ms) 
  b1m = 0.03801123582189022     (/mV) 
  c1m = -3.3336694768115834e-07     (/mV2) 
  d1m = -4.7000105421967445e-09     (/mV3) 
  b2m = 0.03801507790058116     (/mV) 
  c2m = 2.7377623572591333e-07     (/mV2) 
  d2m = -3.339746512763481e-09     (/mV3) 
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
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}