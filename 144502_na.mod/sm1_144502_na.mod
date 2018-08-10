NEURON
{
  SUFFIX iona 
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

  ah = -0.12492522250869244     (/mV) 
  bh = 7.615159061031076     (1) 
  vhh = -64.88933340469359     (mV) 
  Ah = 44.4410504813845     (/ms) 
  b1h = 0.08618552978363303     (/mV) 
  c1h = 0.0007186150862331385     (/mV2) 
  d1h = 9.575519907012326e-06     (/mV3) 
  b2h = 0.05146474016846082     (/mV) 
  c2h = -0.00032400041271215017     (/mV2) 
  d2h = 8.122123183277476e-07     (/mV3) 
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
}

STATE
{
  h
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
}

INITIAL
{
  rates(v)
  h = hInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 


  UNITSON
}