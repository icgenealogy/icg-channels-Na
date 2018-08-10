NEURON
{
  SUFFIX nacurrent 
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

  ah = -0.2500001232687667     (/mV) 
  bh = 12.500007869677631     (1) 
  vhh = -48.23983201336562     (mV) 
  Ah = 18.67185372539621     (/ms) 
  b1h = -0.20935629564590577     (/mV) 
  c1h = 0.002927445600504108     (/mV2) 
  d1h = -1.179012986881505e-05     (/mV3) 
  b2h = -0.2203401893509152     (/mV) 
  c2h = -0.006798187771997567     (/mV2) 
  d2h = -7.333026640573813e-05     (/mV3) 

  aI = -0.005955624414671087     (/mV) 
  bI = -1.6959800980351996     (1) 
  vhI = -65.83432607841861     (mV) 
  AI = 2934.887860352264     (/ms) 
  b1I = -0.37267162605837     (/mV) 
  c1I = 0.001124144397726181     (/mV2) 
  d1I = 7.365128876187846e-06     (/mV3) 
  b2I = -0.08404665702465709     (/mV) 
  c2I = 0.000303424074123373     (/mV2) 
  d2I = 4.832041515029691e-06     (/mV3) 

  am = 0.13888871291147808     (/mV) 
  bm = -5.337840070247862     (1) 
  vhm = -34.00062883133335     (mV) 
  Am = 0.3067083175140181     (/ms) 
  b1m = -0.0789912838024553     (/mV) 
  c1m = 0.0007946009960276957     (/mV2) 
  d1m = -2.5251554718943536e-06     (/mV3) 
  b2m = -0.032045145264008604     (/mV) 
  c2m = 1.6759196426782876e-05     (/mV2) 
  d2m = 2.189165021986104e-06     (/mV3) 
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
  IInf 
  ITau 
  mInf 
  mTau 
}

STATE
{
  h
  I
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  I' = (IInf - I) / ITau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  I = IInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    IInf = 1/(1 + exp(-aI*v + bI)) 
    ITau = AI / ( exp(-(b1I*(v-vhI) + c1I*(v-vhI)^2 + d1I*(v-vhI)^3)) + exp((b2I*(v-vhI) + c2I*(v-vhI)^2 + d2I*(v-vhI)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}