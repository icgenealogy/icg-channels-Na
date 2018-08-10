NEURON
{
  SUFFIX Nasoma 
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

  ah = -0.1414725918767857     (/mV) 
  bh = 8.510933657919312     (1) 
  vhh = -64.40396204903905     (mV) 
  Ah = 17.33712796569069     (/ms) 
  b1h = -0.09377288821068161     (/mV) 
  c1h = 0.0009413613879789784     (/mV2) 
  d1h = -2.944727800950621e-06     (/mV3) 
  b2h = -0.08961057031910775     (/mV) 
  c2h = -0.0019121477154225557     (/mV2) 
  d2h = -2.8966790512887392e-05     (/mV3) 

  am = 0.10469540787755995     (/mV) 
  bm = -3.9338601723688322     (1) 
  vhm = -42.73472631333817     (mV) 
  Am = 0.9705504739503584     (/ms) 
  b1m = -0.0411297560320391     (/mV) 
  c1m = 0.00026418917970783413     (/mV2) 
  d1m = -7.042226876528936e-07     (/mV3) 
  b2m = -0.06915047381999954     (/mV) 
  c2m = -0.0002785567646946672     (/mV2) 
  d2m = -9.36395463803744e-07     (/mV3) 
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