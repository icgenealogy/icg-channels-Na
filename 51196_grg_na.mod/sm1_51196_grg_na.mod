NEURON
{
  SUFFIX GrG_Na 
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

  ah = -0.49868997749602556     (/mV) 
  bh = 18.004973713853364     (1) 
  vhh = -35.0     (mV) 
  Ah = 16.040000001706147     (/ms) 
  b1h = 37.65550851991271     (/mV) 
  c1h = -0.09401831406545211     (/mV2) 
  d1h = 6.396736149994543e-06     (/mV3) 
  b2h = 37.68311188852105     (/mV) 
  c2h = -0.09080959946961269     (/mV2) 
  d2h = -1.5982257866339237e-05     (/mV3) 

  am = 0.1040218176319032     (/mV) 
  bm = -1.9174770299219113     (1) 
  vhm = -28.601657055963496     (mV) 
  Am = 0.06279106991297308     (/ms) 
  b1m = -0.03916189195608382     (/mV) 
  c1m = 0.0004985318004786937     (/mV2) 
  d1m = -2.5286932737276137e-06     (/mV3) 
  b2m = -0.09520006282506475     (/mV) 
  c2m = -0.0029617218470289336     (/mV2) 
  d2m = -4.399581434332897e-05     (/mV3) 
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