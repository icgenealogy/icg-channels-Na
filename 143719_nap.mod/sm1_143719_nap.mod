NEURON
{
  SUFFIX nap 
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

  an = 1.02151808531682     (/mV) 
  bn = -53.08761404789373     (1) 
  vhn = -20.132079849927532     (mV) 
  An = 0.9996749538672128     (/ms) 
  b1n = -0.0015270018210565754     (/mV) 
  c1n = -3.94835000727561e-05     (/mV2) 
  d1n = -2.7787067177699654e-07     (/mV3) 
  b2n = 3.7001067494107907     (/mV) 
  c2n = 0.010820020248510684     (/mV2) 
  d2n = 8.971217606916601e-05     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nInf 
  nTau 
}

STATE
{
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n*n*n
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}