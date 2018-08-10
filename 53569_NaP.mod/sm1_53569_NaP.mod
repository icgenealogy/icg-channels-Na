NEURON
{
  SUFFIX NaP 
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

  ah = -0.14909615866095413     (/mV) 
  bh = 7.444461715044708     (1) 
  vhh = -50.59240086453308     (mV) 
  Ah = 4.814472124946749     (/ms) 
  b1h = -0.10701022713447818     (/mV) 
  c1h = 0.0011560320144803722     (/mV2) 
  d1h = -3.905402536692037e-06     (/mV3) 
  b2h = -0.06925296759473687     (/mV) 
  c2h = -0.0006870426043359371     (/mV2) 
  d2h = -9.04730480129271e-06     (/mV3) 

  am = 0.12544952894987224     (/mV) 
  bm = -3.657437300718694     (1) 
  vhm = -32.10762454015781     (mV) 
  Am = 0.19955120900126735     (/ms) 
  b1m = 0.03081894826515301     (/mV) 
  c1m = -2.8043491151228045e-05     (/mV2) 
  d1m = -2.154517785888116e-06     (/mV3) 
  b2m = 0.03242543494730372     (/mV) 
  c2m = -0.00010956715808708603     (/mV2) 
  d2m = -3.6890719069216505e-07     (/mV3) 
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
  g = gbar*h*m
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