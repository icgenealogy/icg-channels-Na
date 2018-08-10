NEURON
{
  SUFFIX na_sej 
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

  ah = -0.16129015130014016     (/mV) 
  bh = 10.483866196061618     (1) 
  vhh = -61.839875919644555     (mV) 
  Ah = 23.3612661230549     (/ms) 
  b1h = 0.1544151666452644     (/mV) 
  c1h = 0.002974543782633336     (/mV2) 
  d1h = 2.6469428450243305e-05     (/mV3) 
  b2h = 0.178624247680912     (/mV) 
  c2h = -0.002291580592229332     (/mV2) 
  d2h = 9.905893255581303e-06     (/mV3) 

  am = 0.11111084561240928     (/mV) 
  bm = -4.272599623781179     (1) 
  vhm = -37.40047439756125     (mV) 
  Am = 0.18059965188579283     (/ms) 
  b1m = -0.0451935487165594     (/mV) 
  c1m = 0.00020198382306724872     (/mV2) 
  d1m = 3.640698054055314e-07     (/mV3) 
  b2m = -0.0349481194704929     (/mV) 
  c2m = 3.148663488342266e-05     (/mV2) 
  d2m = 2.6799513201909388e-06     (/mV3) 
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