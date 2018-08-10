NEURON
{
  SUFFIX Golgi_NaR 
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

  as = 0.03439689606679359     (/mV) 
  bs = 0.7669055627121748     (1) 
  vhs = -46.34091693834107     (mV) 
  As = 4.132416230594445     (/ms) 
  b1s = 0.005364344251066176     (/mV) 
  c1s = -0.0006231018629830464     (/mV2) 
  d1s = 4.208828911892106e-06     (/mV3) 
  b2s = -0.1417603689202406     (/mV) 
  c2s = 0.002054434522529298     (/mV2) 
  d2s = -6.960926624686603e-06     (/mV3) 

  af = -0.07825769977237026     (/mV) 
  bf = 3.0225648745148592     (1) 
  vhf = -33.059133457198094     (mV) 
  Af = 0.8292020464793795     (/ms) 
  b1f = 0.005538883804344022     (/mV) 
  c1f = -0.0001520043532014586     (/mV2) 
  d1f = -7.724163737987582e-07     (/mV3) 
  b2f = 0.06727896515918723     (/mV) 
  c2f = -1.5593983213712613e-05     (/mV2) 
  d2f = -2.8466772574096397e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  sInf 
  sTau 
  fInf 
  fTau 
}

STATE
{
  s
  f
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*s*f
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  s' = (sInf - s) / sTau 
  f' = (fInf - f) / fTau 
}

INITIAL
{
  rates(v)
  s = sInf 
  f = fInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    sInf = 1/(1 + exp(-as*v + bs)) 
    sTau = As / ( exp(-(b1s*(v-vhs) + c1s*(v-vhs)^2 + d1s*(v-vhs)^3)) + exp((b2s*(v-vhs) + c2s*(v-vhs)^2 + d2s*(v-vhs)^3)) ) 

    fInf = 1/(1 + exp(-af*v + bf)) 
    fTau = Af / ( exp(-(b1f*(v-vhf) + c1f*(v-vhf)^2 + d1f*(v-vhf)^3)) + exp((b2f*(v-vhf) + c2f*(v-vhf)^2 + d2f*(v-vhf)^3)) ) 


  UNITSON
}