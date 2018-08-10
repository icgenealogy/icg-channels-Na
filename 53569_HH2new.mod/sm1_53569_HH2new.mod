NEURON
{
  SUFFIX HH2new 
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

  am = 0.11295691107730477     (/mV) 
  bm = -3.3143106862649745     (1) 
  vhm = -55.53154553611752     (mV) 
  Am = 0.14642463860939908     (/ms) 
  b1m = -0.008824822922760848     (/mV) 
  c1m = -0.0001589437092559967     (/mV2) 
  d1m = 1.0556421385405217e-06     (/mV3) 
  b2m = -0.0375419775899118     (/mV) 
  c2m = -0.00020494410557330244     (/mV2) 
  d2m = 3.2672772603249482e-06     (/mV3) 

  ah = -0.2474886938383068     (/mV) 
  bh = 8.246265533493924     (1) 
  vhh = -33.65062596885506     (mV) 
  Ah = 9.095558867761355     (/ms) 
  b1h = -0.20232871231752697     (/mV) 
  c1h = 0.002929881382337901     (/mV2) 
  d1h = -1.2501029411081771e-05     (/mV3) 
  b2h = -0.0612034292255706     (/mV) 
  c2h = -0.00026210402628733655     (/mV2) 
  d2h = -2.9900381116103855e-06     (/mV3) 

  an = 0.09068689383249524     (/mV) 
  bn = -2.611512372187617     (1) 
  vhn = -51.29893418383732     (mV) 
  An = 2.886790249536415     (/ms) 
  b1n = -0.04191407570642935     (/mV) 
  c1n = 0.0002739553551948287     (/mV2) 
  d1n = -7.539535344990324e-07     (/mV3) 
  b2n = -0.06990263673383365     (/mV) 
  c2n = -0.0011921993227104525     (/mV2) 
  d2n = -1.096148589889475e-05     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
  hInf 
  hTau 
  nInf 
  nTau 
}

STATE
{
  m
  h
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*m*h*n*n*n*n
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  m = mInf 
  h = hInf 
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}