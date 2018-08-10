NEURON
{
  SUFFIX NaxSH10_ChannelML 
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

  ah = -0.250001130651872     (/mV) 
  bh = 10.000050222284282     (1) 
  vhh = -37.99366210903517     (mV) 
  Ah = 14.990354456463264     (/ms) 
  b1h = 0.203255753664351     (/mV) 
  c1h = 0.0054701046730137115     (/mV2) 
  d1h = 5.0749035462151924e-05     (/mV3) 
  b2h = 0.21196855018618757     (/mV) 
  c2h = -0.003170634805510527     (/mV2) 
  d2h = 1.3631211766045518e-05     (/mV3) 

  am = 0.13888850827297494     (/mV) 
  bm = -3.9489512416830967     (1) 
  vhm = -21.159046245834187     (mV) 
  Am = 0.23231397461201367     (/ms) 
  b1m = -0.08416249984370827     (/mV) 
  c1m = 0.0009623368655569791     (/mV2) 
  d1m = -2.9110814763290898e-06     (/mV3) 
  b2m = -0.023975471219545366     (/mV) 
  c2m = 0.00010476621845080575     (/mV2) 
  d2m = 2.0235698017719897e-06     (/mV3) 
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