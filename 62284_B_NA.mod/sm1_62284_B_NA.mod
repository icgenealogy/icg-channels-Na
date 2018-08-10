NEURON
{
  SUFFIX B_Na 
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

  ah = -0.11111120763956739     (/mV) 
  bh = 7.111118879369787     (1) 
  vhh = -43.69560367378991     (mV) 
  Ah = 1.116593137752104     (/ms) 
  b1h = 0.015865438702680576     (/mV) 
  c1h = -4.244351228039226e-05     (/mV2) 
  d1h = -2.735210060761516e-07     (/mV3) 
  b2h = 0.23394921029543567     (/mV) 
  c2h = -0.0045288557731630015     (/mV2) 
  d2h = 2.9205473616803632e-05     (/mV3) 

  am = 0.11111174288314962     (/mV) 
  bm = -3.4948864714296057     (1) 
  vhm = -43.38569340179288     (mV) 
  Am = 1.0914715540851774     (/ms) 
  b1m = 0.01470409178615365     (/mV) 
  c1m = -2.0317533357820558e-05     (/mV2) 
  d1m = 5.555189683969719e-07     (/mV3) 
  b2m = 0.24843988310654203     (/mV) 
  c2m = -0.005888129131625044     (/mV2) 
  d2m = 5.307842226736303e-05     (/mV3) 
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