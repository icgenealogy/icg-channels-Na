NEURON
{
  SUFFIX na 
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

  ah = -0.1666678610592848     (/mV) 
  bh = 9.333430614939555     (1) 
  vhh = -31.384270884872315     (mV) 
  Ah = 19.642926796751038     (/ms) 
  b1h = -0.10171171295137121     (/mV) 
  c1h = -0.0005672039269962811     (/mV2) 
  d1h = 1.7293264727206175e-05     (/mV3) 
  b2h = 0.0184380186254407     (/mV) 
  c2h = 0.0007544709848511775     (/mV2) 
  d2h = -5.051249364838695e-06     (/mV3) 

  am = 0.09999954765035772     (/mV) 
  bm = -3.3000083035672674     (1) 
  vhm = -26.029748068175536     (mV) 
  Am = 0.3501233419864821     (/ms) 
  b1m = -0.06162721149032775     (/mV) 
  c1m = 0.0006937790499144983     (/mV2) 
  d1m = -2.5219353829320058e-06     (/mV3) 
  b2m = -0.04314983135788186     (/mV) 
  c2m = -0.00017660500619791472     (/mV2) 
  d2m = 1.1563930146679165e-06     (/mV3) 
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