NEURON
{
  SUFFIX na12 
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

  ah = -0.16110853117691531     (/mV) 
  bh = 10.7942638544828     (1) 
  vhh = -61.80136064804162     (mV) 
  Ah = 19.212436790059513     (/ms) 
  b1h = 0.14671921383930492     (/mV) 
  c1h = 0.0029002287811235937     (/mV2) 
  d1h = 2.6309971770612447e-05     (/mV3) 
  b2h = 0.14185716626520772     (/mV) 
  c2h = -0.0014991413898000196     (/mV2) 
  d2h = 5.4164800880267446e-06     (/mV3) 

  am = 0.14273535732448447     (/mV) 
  bm = -5.807382291595278     (1) 
  vhm = -45.79142554012623     (mV) 
  Am = 0.2977622451254423     (/ms) 
  b1m = -0.05614518226414453     (/mV) 
  c1m = 0.00041742780829248687     (/mV2) 
  d1m = -1.1158730244024082e-06     (/mV3) 
  b2m = -0.07514624526427317     (/mV) 
  c2m = -0.0010534977443832914     (/mV2) 
  d2m = -5.231890566595067e-06     (/mV3) 
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