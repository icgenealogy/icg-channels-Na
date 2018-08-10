NEURON
{
  SUFFIX nainter 
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

  ah = -0.11226945803715291     (/mV) 
  bh = 3.25582477186598     (1) 
  vhh = -29.967456402408793     (mV) 
  Ah = 22.764588722897678     (/ms) 
  b1h = -0.07829018437509329     (/mV) 
  c1h = -9.345285715211456e-05     (/mV2) 
  d1h = 9.015976828761367e-07     (/mV3) 
  b2h = -0.03280294451602887     (/mV) 
  c2h = -4.297165853117053e-05     (/mV2) 
  d2h = -2.537488972137963e-07     (/mV3) 

  am = 0.12214376531418508     (/mV) 
  bm = -2.2596633062852955     (1) 
  vhm = 21.669258093858794     (mV) 
  Am = 0.28548288158289514     (/ms) 
  b1m = 0.23986697899085896     (/mV) 
  c1m = 0.009615004658669829     (/mV2) 
  d1m = 6.589649716613464e-05     (/mV3) 
  b2m = 0.037671970466583465     (/mV) 
  c2m = 0.0010278669429396405     (/mV2) 
  d2m = 5.025991164240495e-06     (/mV3) 
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
  g = gbar*h*m*m*m
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