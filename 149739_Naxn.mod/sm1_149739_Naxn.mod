NEURON
{
  SUFFIX nax 
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

  ah = -0.24999986821222347     (/mV) 
  bh = 11.249990779587572     (1) 
  vhh = -43.293511283648186     (mV) 
  Ah = 15.019292396155768     (/ms) 
  b1h = -0.20501905702887413     (/mV) 
  c1h = 0.0029560490907511996     (/mV2) 
  d1h = -1.2243667529773793e-05     (/mV3) 
  b2h = -0.21404785663394552     (/mV) 
  c2h = -0.006233303413335527     (/mV2) 
  d2h = -6.27624613291557e-05     (/mV3) 

  am = 0.13888872874648067     (/mV) 
  bm = -4.643393337872212     (1) 
  vhm = -26.23369014615548     (mV) 
  Am = 0.23275517095574202     (/ms) 
  b1m = 0.02495525081638534     (/mV) 
  c1m = -6.85178526120918e-05     (/mV2) 
  d1m = -1.6777384085493585e-06     (/mV3) 
  b2m = 0.08535469101772408     (/mV) 
  c2m = -0.001028955286140716     (/mV2) 
  d2m = 3.8144646926022276e-06     (/mV3) 
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