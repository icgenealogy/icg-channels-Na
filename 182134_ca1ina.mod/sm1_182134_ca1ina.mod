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

  ai = -0.005955739458234923     (/mV) 
  bi = -1.6959838242727767     (1) 
  vhi = -65.83432607841861     (mV) 
  Ai = 2934.887860352264     (/ms) 
  b1i = -0.37267162605837     (/mV) 
  c1i = 0.001124144397726181     (/mV2) 
  d1i = 7.365128876187846e-06     (/mV3) 
  b2i = -0.08404665702465709     (/mV) 
  c2i = 0.000303424074123373     (/mV2) 
  d2i = 4.832041515029691e-06     (/mV3) 

  ah = -0.2500001232687576     (/mV) 
  bh = 12.500007869677004     (1) 
  vhh = -48.23983201336562     (mV) 
  Ah = 18.67185372539621     (/ms) 
  b1h = -0.20935629564590577     (/mV) 
  c1h = 0.002927445600504108     (/mV2) 
  d1h = -1.179012986881505e-05     (/mV3) 
  b2h = -0.2203401893509152     (/mV) 
  c2h = -0.006798187771997567     (/mV2) 
  d2h = -7.333026640573813e-05     (/mV3) 

  am = 0.13888871291144977     (/mV) 
  bm = -5.337840070246734     (1) 
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
  iInf 
  iTau 
  hInf 
  hTau 
  mInf 
  mTau 
}

STATE
{
  i
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*i*h*m*m*m
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  i' = (iInf - i) / iTau 
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  i = iInf 
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    iInf = 1/(1 + exp(-ai*v + bi)) 
    iTau = Ai / ( exp(-(b1i*(v-vhi) + c1i*(v-vhi)^2 + d1i*(v-vhi)^3)) + exp((b2i*(v-vhi) + c2i*(v-vhi)^2 + d2i*(v-vhi)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}