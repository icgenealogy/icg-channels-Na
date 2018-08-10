NEURON
{
  SUFFIX B_Na 
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

  ah = -0.11111096804638848     (/mV) 
  bh = 8.333323855439549     (1) 
  vhh = -67.42542176107312     (mV) 
  Ah = 2.1169279614689986     (/ms) 
  b1h = 0.1579444059852811     (/mV) 
  c1h = 0.0037793782931796146     (/mV2) 
  d1h = 4.4290304045573824e-05     (/mV3) 
  b2h = 0.1345905160914471     (/mV) 
  c2h = -0.001349070574731392     (/mV2) 
  d2h = 4.468578050987078e-06     (/mV3) 

  am = 0.11111097389862959     (/mV) 
  bm = -4.272607197672756     (1) 
  vhm = -65.83706616362437     (mV) 
  Am = 2.121104938171783     (/ms) 
  b1m = 0.12765445327839825     (/mV) 
  c1m = 0.0017381865265018765     (/mV2) 
  d1m = 7.528094289446758e-06     (/mV3) 
  b2m = 0.14247344725148792     (/mV) 
  c2m = -0.0017930288820020722     (/mV2) 
  d2m = 9.011978197189056e-06     (/mV3) 
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