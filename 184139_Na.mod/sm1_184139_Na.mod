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

  ah = -0.12499996976283675     (/mV) 
  bh = 9.125008249306237     (1) 
  vhh = -66.92423838229683     (mV) 
  Ah = 34.47504413159874     (/ms) 
  b1h = -0.14079609836269041     (/mV) 
  c1h = 0.0014629776356923033     (/mV2) 
  d1h = -5.1794635212082924e-06     (/mV3) 
  b2h = -0.1499480575298776     (/mV) 
  c2h = -0.0031514599041875366     (/mV2) 
  d2h = -3.1810164393399415e-05     (/mV3) 

  am = 0.12499983710501615     (/mV) 
  bm = -5.1337144970877695     (1) 
  vhm = -47.63801221050546     (mV) 
  Am = 0.46389227607173905     (/ms) 
  b1m = -0.048227469586540624     (/mV) 
  c1m = 0.00030284941594914324     (/mV2) 
  d1m = -7.08954530032178e-07     (/mV3) 
  b2m = -0.0712432430677494     (/mV) 
  c2m = -0.0011085848909195867     (/mV2) 
  d2m = -7.69037064914426e-06     (/mV3) 
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