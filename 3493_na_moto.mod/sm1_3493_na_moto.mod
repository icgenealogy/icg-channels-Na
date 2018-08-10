NEURON
{
  SUFFIX namot 
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

  ah = -0.09803912795655022     (/mV) 
  bh = 7.999996579385781     (1) 
  vhh = -2.084830544263585     (mV) 
  Ah = 0.4495590415585725     (/ms) 
  b1h = -0.02752885323436011     (/mV) 
  c1h = 0.0003653044796870123     (/mV2) 
  d1h = -1.5726442572132167e-06     (/mV3) 
  b2h = -0.007592224205930244     (/mV) 
  c2h = -0.0011570947283239942     (/mV2) 
  d2h = -9.98105370103097e-06     (/mV3) 

  am = 0.16393415437177764     (/mV) 
  bm = -6.377034570340227     (1) 
  vhm = -32.60094442945868     (mV) 
  Am = 0.3789879614912129     (/ms) 
  b1m = 0.01687900064256622     (/mV) 
  c1m = 0.00028426381643850446     (/mV2) 
  d1m = 2.620812381368759e-06     (/mV3) 
  b2m = 0.05777746817396129     (/mV) 
  c2m = -0.00031419613134276055     (/mV2) 
  d2m = 1.6282294911997766e-06     (/mV3) 
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