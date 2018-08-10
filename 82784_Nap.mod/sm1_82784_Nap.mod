NEURON
{
  SUFFIX Nap 
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

  ah = -0.2549304736386867     (/mV) 
  bh = 14.369960628091635     (1) 
  vhh = -56.364185685076926     (mV) 
  Ah = 1242.7775187270563     (/ms) 
  b1h = -0.006467526092021201     (/mV) 
  c1h = 7.507730585572684e-07     (/mV2) 
  d1h = 2.8565957338665446e-09     (/mV3) 
  b2h = -0.2484801221019655     (/mV) 
  c2h = -2.7567398882985275e-06     (/mV2) 
  d2h = 1.3000227084572642e-08     (/mV3) 

  am = 0.09667127739120644     (/mV) 
  bm = -0.2006496408163099     (1) 
  vhm = 27.669942969120818     (mV) 
  Am = 0.18892610041830857     (/ms) 
  b1m = -0.039058402354734144     (/mV) 
  c1m = 0.0005523869927413765     (/mV2) 
  d1m = -3.0133567521047285e-06     (/mV3) 
  b2m = -0.006331834198615209     (/mV) 
  c2m = 0.00015372917054454256     (/mV2) 
  d2m = 9.440286915224438e-07     (/mV3) 
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