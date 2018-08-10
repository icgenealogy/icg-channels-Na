NEURON
{
  SUFFIX NaTa_t 
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

  ah = -0.16666639513868026     (/mV) 
  bh = 10.999981729139636     (1) 
  vhh = -72.16227780068864     (mV) 
  Ah = 3.483070457202446     (/ms) 
  b1h = -0.05460886998401727     (/mV) 
  c1h = 0.0003834132386300533     (/mV2) 
  d1h = -1.0306008963353065e-06     (/mV3) 
  b2h = -0.10671579781129222     (/mV) 
  c2h = -0.0027020525903908264     (/mV2) 
  d2h = -3.538664666566817e-05     (/mV3) 

  am = 0.1666663306077527     (/mV) 
  bm = -6.717038789860857     (1) 
  vhm = -44.85863709900642     (mV) 
  Am = 0.37947077237183335     (/ms) 
  b1m = -0.06677529926511423     (/mV) 
  c1m = 0.0005880171190539407     (/mV2) 
  d1m = -1.9005375667878346e-06     (/mV3) 
  b2m = -0.09395067064249205     (/mV) 
  c2m = -0.001831374078141017     (/mV2) 
  d2m = -1.4751415073302824e-05     (/mV3) 
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