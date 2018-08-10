NEURON
{
  SUFFIX NaFax 
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

  ah = -0.2340551504443316     (/mV) 
  bh = 8.439337201860138     (1) 
  vhh = -41.88337867379786     (mV) 
  Ah = 1.781077324349133     (/ms) 
  b1h = 0.12296681371581968     (/mV) 
  c1h = 0.0027350282436633703     (/mV2) 
  d1h = 3.188400909352168e-05     (/mV3) 
  b2h = 0.13665070725662556     (/mV) 
  c2h = -0.001791607223753524     (/mV2) 
  d2h = 7.015111127386209e-06     (/mV3) 

  am = 0.13999698771162197     (/mV) 
  bm = -5.00688881707605     (1) 
  vhm = -35.323549387193346     (mV) 
  Am = 0.10639463534696995     (/ms) 
  b1m = 0.04225594177781041     (/mV) 
  c1m = 0.00048317726583443573     (/mV2) 
  d1m = 3.4322837761320306e-06     (/mV3) 
  b2m = 0.05098044714710637     (/mV) 
  c2m = -0.000500944231722518     (/mV2) 
  d2m = 1.5956986788949033e-06     (/mV3) 
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