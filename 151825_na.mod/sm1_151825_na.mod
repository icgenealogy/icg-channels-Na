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

  ah = -0.16666638818679602     (/mV) 
  bh = 9.333315962527408     (1) 
  vhh = -20.722032147311154     (mV) 
  Ah = 2.3055351264025035     (/ms) 
  b1h = -0.13002418996759008     (/mV) 
  c1h = 0.0015108072008815816     (/mV2) 
  d1h = -5.442274695844736e-06     (/mV3) 
  b2h = 0.04239617895193789     (/mV) 
  c2h = 0.0004300417021508381     (/mV2) 
  d2h = -8.224637907204083e-06     (/mV3) 

  am = 0.09999981671158166     (/mV) 
  bm = -3.2999885526798955     (1) 
  vhm = -25.31593650631692     (mV) 
  Am = 0.10126681417763557     (/ms) 
  b1m = 0.04276055499353669     (/mV) 
  c1m = 0.00012896204539781432     (/mV2) 
  d1m = -2.0966876026142813e-06     (/mV3) 
  b2m = 0.06271253848959868     (/mV) 
  c2m = -0.0007623142450167952     (/mV2) 
  d2m = 2.936574287257792e-06     (/mV3) 
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