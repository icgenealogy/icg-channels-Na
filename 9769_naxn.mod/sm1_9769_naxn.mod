NEURON
{
  SUFFIX nax 
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

  ah = -0.24999955331498572     (/mV) 
  bh = 8.749977903508608     (1) 
  vhh = -32.699781887608694     (mV) 
  Ah = 14.93356861207547     (/ms) 
  b1h = 0.19341853370357862     (/mV) 
  c1h = 0.004834030418191378     (/mV2) 
  d1h = 4.154002411356305e-05     (/mV3) 
  b2h = 0.21900498226065818     (/mV) 
  c2h = -0.003403357442411117     (/mV2) 
  d2h = 1.5213382296675904e-05     (/mV3) 

  am = 0.1388886543112195     (/mV) 
  bm = -3.2545026218941766     (1) 
  vhm = -19.84008685606428     (mV) 
  Am = 0.2559949006588593     (/ms) 
  b1m = -0.08524014108112386     (/mV) 
  c1m = 0.0009936967325043196     (/mV2) 
  d1m = -3.6346429087031654e-06     (/mV3) 
  b2m = -0.040052154602597534     (/mV) 
  c2m = -0.0003153325479903794     (/mV2) 
  d2m = -1.1527374857194237e-06     (/mV3) 
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