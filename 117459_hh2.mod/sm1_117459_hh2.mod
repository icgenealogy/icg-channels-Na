NEURON
{
  SUFFIX hh2 
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

  ah = -0.11530037460467249     (/mV) 
  bh = 6.650258796535576     (1) 
  vhh = -64.9575794428315     (mV) 
  Ah = 0.6481557604622592     (/ms) 
  b1h = -0.08986714912241928     (/mV) 
  c1h = 0.0008862100948655642     (/mV2) 
  d1h = -2.7286280410374067e-06     (/mV3) 
  b2h = -0.0660864732080897     (/mV) 
  c2h = -0.0016466479927610715     (/mV2) 
  d2h = -2.199417718556478e-05     (/mV3) 

  an = 0.05115431598302338     (/mV) 
  bn = -2.6801853090300463     (1) 
  vhn = -74.83658172249544     (mV) 
  An = 0.3777939565415886     (/ms) 
  b1n = -0.008394430247404182     (/mV) 
  c1n = -0.0003096819696467863     (/mV2) 
  d1n = 3.074916534099031e-06     (/mV3) 
  b2n = -1.0852075476035141e-09     (/mV) 
  c2n = 0.0002951842793028292     (/mV2) 
  d2n = -1.219421729339401e-06     (/mV3) 

  am = 0.09830793111954274     (/mV) 
  bm = -3.57316352456176     (1) 
  vhm = -67.20022630522016     (mV) 
  Am = 0.6548943079088829     (/ms) 
  b1m = 0.08117328731929696     (/mV) 
  c1m = 0.0032661759428040567     (/mV2) 
  d1m = 6.791466910729105e-05     (/mV3) 
  b2m = 0.08364381188297675     (/mV) 
  c2m = -0.0007544761876659573     (/mV2) 
  d2m = 2.1359226057233914e-06     (/mV3) 
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
  nInf 
  nTau 
  mInf 
  mTau 
}

STATE
{
  h
  n
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*n*n*n*n*m*m*m
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  n = nInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}