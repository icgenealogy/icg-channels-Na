NEURON
{
  SUFFIX nahh 
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

  ah = -0.2475033923801006     (/mV) 
  bh = 10.721938738633193     (1) 
  vhh = -44.33305037035574     (mV) 
  Ah = 4.915571094515876     (/ms) 
  b1h = -0.19502846365281137     (/mV) 
  c1h = 0.002648292355057956     (/mV2) 
  d1h = -1.0545995808948608e-05     (/mV3) 
  b2h = -0.07174629015519635     (/mV) 
  c2h = -0.0006856235091814014     (/mV2) 
  d2h = -8.003226865106709e-06     (/mV3) 

  am = 0.13491621162181613     (/mV) 
  bm = -5.251774156071185     (1) 
  vhm = -52.58260101386403     (mV) 
  Am = 0.11277997372225483     (/ms) 
  b1m = -0.024495147095919757     (/mV) 
  c1m = 2.0890054456916528e-06     (/mV2) 
  d1m = 7.064203502368652e-07     (/mV3) 
  b2m = -0.043771197559281216     (/mV) 
  c2m = -0.00019111780122931936     (/mV2) 
  d2m = 3.790277029145558e-06     (/mV3) 
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