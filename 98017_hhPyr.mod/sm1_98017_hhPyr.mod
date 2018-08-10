NEURON
{
  SUFFIX hhPyr 
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

  ah = -0.1414767959860261     (/mV) 
  bh = 6.2475503500652865     (1) 
  vhh = -46.44780663899414     (mV) 
  Ah = 0.5956064370458408     (/ms) 
  b1h = -0.09661633951680482     (/mV) 
  c1h = 0.0010626022295823457     (/mV2) 
  d1h = -3.659565732338684e-06     (/mV3) 
  b2h = -0.06860034448189883     (/mV) 
  c2h = -0.000708505020809205     (/mV2) 
  d2h = -7.311275468750544e-06     (/mV3) 

  an = 0.05720162611790115     (/mV) 
  bn = -1.71522766922599     (1) 
  vhn = -52.15867234439556     (mV) 
  An = 0.3987513599935466     (/ms) 
  b1n = -0.04119791493301328     (/mV) 
  c1n = 0.0002637015363414731     (/mV2) 
  d1n = -7.031312958568088e-07     (/mV3) 
  b2n = -0.03537453709688636     (/mV) 
  c2n = -0.0004997733369558299     (/mV2) 
  d2n = -4.308912534152548e-06     (/mV3) 
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
}

STATE
{
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*n*n*n*n
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}