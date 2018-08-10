NEURON
{
  SUFFIX napf 
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

  am = 0.10018028852724489     (/mV) 
  bm = -3.8092800821870876     (1) 
  vhm = -33.068778479253865     (mV) 
  Am = 0.28515275972910026     (/ms) 
  b1m = 0.1767946201224783     (/mV) 
  c1m = 0.004436021252694326     (/mV2) 
  d1m = 3.5045107615474516e-05     (/mV3) 
  b2m = 0.10821890952526844     (/mV) 
  c2m = -0.001477520422015668     (/mV2) 
  d2m = 6.076812760223376e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*m
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}