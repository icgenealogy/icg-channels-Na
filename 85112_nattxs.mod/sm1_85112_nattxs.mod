NEURON
{
  SUFFIX nattxs 
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

  ah = -0.104169336537098     (/mV) 
  bh = 5.271480277287399     (1) 
  vhh = -46.85672937713155     (mV) 
  Ah = 14.830954335131661     (/ms) 
  b1h = -0.09285299726701121     (/mV) 
  c1h = 0.00021136356549974252     (/mV2) 
  d1h = 1.9029643873339443e-06     (/mV3) 
  b2h = -0.010900738217184723     (/mV) 
  c2h = 0.0001788870634592608     (/mV2) 
  d2h = 9.816457335786682e-07     (/mV3) 

  am = 0.11524857456157923     (/mV) 
  bm = -3.3211407896086147     (1) 
  vhm = -35.556464286865264     (mV) 
  Am = 0.4633381660680909     (/ms) 
  b1m = 0.07184500427965651     (/mV) 
  c1m = 0.0011358092009148715     (/mV2) 
  d1m = 7.343983130120015e-06     (/mV3) 
  b2m = 0.06839482387924284     (/mV) 
  c2m = -0.0007365662152242756     (/mV2) 
  d2m = 2.554834341197245e-06     (/mV3) 
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