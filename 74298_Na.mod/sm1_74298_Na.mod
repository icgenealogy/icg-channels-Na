NEURON
{
  SUFFIX Na 
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

  ah = -0.24748573958540124     (/mV) 
  bh = 9.483506649887092     (1) 
  vhh = -38.986897954585025     (mV) 
  Ah = 6.391689546472675     (/ms) 
  b1h = 0.06642008530736916     (/mV) 
  c1h = 0.00046492742828004355     (/mV2) 
  d1h = 5.296154001601246e-06     (/mV3) 
  b2h = 0.1986275657947956     (/mV) 
  c2h = -0.0027859880887974274     (/mV2) 
  d2h = 1.1481202149654018e-05     (/mV3) 

  am = 0.16139095640555665     (/mV) 
  bm = -5.65586011242585     (1) 
  vhm = -72.5260780150865     (mV) 
  Am = 0.09013261122138436     (/ms) 
  b1m = 0.017572631535586674     (/mV) 
  c1m = 0.0013963599973162814     (/mV2) 
  d1m = -1.0479634982802834e-05     (/mV3) 
  b2m = -0.01756996904356176     (/mV) 
  c2m = 0.000545188759864853     (/mV2) 
  d2m = -2.6040330433012236e-06     (/mV3) 
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
  g = gbar*h*m*m
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