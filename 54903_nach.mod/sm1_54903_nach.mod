NEURON
{
  SUFFIX NaCh 
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

  ah = -0.15999985332532857     (/mV) 
  bh = 9.640614341343909     (1) 
  vhh = -60.167236712795685     (mV) 
  Ah = 5.085544322798033     (/ms) 
  b1h = -0.03839615216965504     (/mV) 
  c1h = -5.460437677400802e-06     (/mV2) 
  d1h = 1.0452078227901798e-07     (/mV3) 
  b2h = -0.12035252051338508     (/mV) 
  c2h = 0.00011645673059300878     (/mV2) 
  d2h = 4.2495863277704205e-06     (/mV3) 

  am = 0.07999981741956322     (/mV) 
  bm = -2.3999911428300567     (1) 
  vhm = -60.0     (mV) 
  Am = 6.0200000000000005     (/ms) 
  b1m = 121.35510706724152     (/mV) 
  c1m = 0.010256764167028289     (/mV2) 
  d1m = 4.889242262963112e-06     (/mV3) 
  b2m = 121.3780770928021     (/mV) 
  c2m = 0.011785404982679937     (/mV2) 
  d2m = -4.084533917234969e-06     (/mV3) 
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