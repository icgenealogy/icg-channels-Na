NEURON
{
  SUFFIX sdsr 
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

  aasd = 0.09999993775348366     (/mV) 
  basd = -3.999996921666495     (1) 
  vhasd = -46.69413497772626     (mV) 
  Aasd = 5.360000002236288     (/ms) 
  b1asd = 1.7694449145534415e-08     (/mV) 
  c1asd = 4.070021865019931e-09     (/mV2) 
  d1asd = 3.239967056968791e-11     (/mV3) 
  b2asd = 1.7734839958276396e-08     (/mV) 
  c2asd = 4.070860013811846e-09     (/mV2) 
  d2asd = 3.2381573228246595e-11     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  asdInf 
  asdTau 
}

STATE
{
  asd
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*asd
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  asd' = (asdInf - asd) / asdTau 
}

INITIAL
{
  rates(v)
  asd = asdInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    asdInf = 1/(1 + exp(-aasd*v + basd)) 
    asdTau = Aasd / ( exp(-(b1asd*(v-vhasd) + c1asd*(v-vhasd)^2 + d1asd*(v-vhasd)^3)) + exp((b2asd*(v-vhasd) + c2asd*(v-vhasd)^2 + d2asd*(v-vhasd)^3)) ) 


  UNITSON
}