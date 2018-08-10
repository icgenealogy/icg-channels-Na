NEURON
{
  SUFFIX CA1NaP 
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

  amP = 0.1999999254383508     (/mV) 
  bmP = -9.799993183208587     (1) 
  vhmP = -18.037771471201047     (mV) 
  AmP = 0.2778915009415232     (/ms) 
  b1mP = -0.06853078789281809     (/mV) 
  c1mP = 0.0007155444545204153     (/mV2) 
  d1mP = -4.077914506927146e-06     (/mV3) 
  b2mP = -0.0018201186136783864     (/mV) 
  c2mP = 7.638281832090636e-05     (/mV2) 
  d2mP = 2.134523896340891e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mPInf 
  mPTau 
}

STATE
{
  mP
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*mP
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  mP' = (mPInf - mP) / mPTau 
}

INITIAL
{
  rates(v)
  mP = mPInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mPInf = 1/(1 + exp(-amP*v + bmP)) 
    mPTau = AmP / ( exp(-(b1mP*(v-vhmP) + c1mP*(v-vhmP)^2 + d1mP*(v-vhmP)^3)) + exp((b2mP*(v-vhmP) + c2mP*(v-vhmP)^2 + d2mP*(v-vhmP)^3)) ) 


  UNITSON
}