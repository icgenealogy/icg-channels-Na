NEURON
{
  SUFFIX wb 
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

  ah = -0.14808174753132505     (/mV) 
  bh = 10.507664331281465     (1) 
  vhh = -64.44536176572437     (mV) 
  Ah = 57.1488819940718     (/ms) 
  b1h = 0.015228189902133165     (/mV) 
  c1h = -0.0004731276890236257     (/mV2) 
  d1h = 1.6135860031138372e-06     (/mV3) 
  b2h = 0.13138109860802263     (/mV) 
  c2h = -0.00131168585532646     (/mV2) 
  d2h = 3.980999282911704e-06     (/mV3) 

  an = 0.05720161673229479     (/mV) 
  bn = -1.7152252288559764     (1) 
  vhn = -55.78613032564703     (mV) 
  An = 11.603096372276267     (/ms) 
  b1n = 0.037907988762120606     (/mV) 
  c1n = 0.000490484436587114     (/mV2) 
  d1n = 3.6835678852603233e-06     (/mV3) 
  b2n = 0.03842569056539727     (/mV) 
  c2n = -0.00020458289242483976     (/mV2) 
  d2n = 4.592941678599075e-07     (/mV3) 
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