NEURON
{
  SUFFIX Naf 
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

  ah = -0.1414762491654464     (/mV) 
  bh = 5.6957670235647     (1) 
  vhh = -41.35768190813034     (mV) 
  Ah = 11.923114466289523     (/ms) 
  b1h = 0.06359845484797201     (/mV) 
  c1h = 0.0005621853109782146     (/mV2) 
  d1h = 6.280000442424234e-06     (/mV3) 
  b2h = 0.10329514384625549     (/mV) 
  c2h = -0.0011746098510093054     (/mV2) 
  d2h = 4.19108865745068e-06     (/mV3) 

  am = 0.0966714070151006     (/mV) 
  bm = -1.7473993237951415     (1) 
  vhm = -29.47526489831788     (mV) 
  Am = 0.2243602552703792     (/ms) 
  b1m = 0.031965904374381834     (/mV) 
  c1m = 5.8527183893779365e-05     (/mV2) 
  d1m = -3.069187202888886e-06     (/mV3) 
  b2m = 0.013299791036635608     (/mV) 
  c2m = 0.00021221312256497534     (/mV2) 
  d2m = -1.8718728365240889e-06     (/mV3) 
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