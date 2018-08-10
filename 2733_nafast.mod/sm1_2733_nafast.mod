NEURON
{
  SUFFIX nafast 
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

  ah = -0.2474909062025882     (/mV) 
  bh = 8.246333591617791     (1) 
  vhh = -33.68144965093311     (mV) 
  Ah = 10.174641470488247     (/ms) 
  b1h = -0.20228892049398645     (/mV) 
  c1h = 0.0029230897429134333     (/mV2) 
  d1h = -1.2464696321975868e-05     (/mV3) 
  b2h = -0.061783311673527214     (/mV) 
  c2h = -0.0002857940963948566     (/mV2) 
  d2h = -3.2637739522152687e-06     (/mV3) 

  am = 0.13491561287719586     (/mV) 
  bm = -3.9160099128275005     (1) 
  vhm = -43.60619070311865     (mV) 
  Am = 0.223808998645513     (/ms) 
  b1m = 0.05629624624786105     (/mV) 
  c1m = 0.0007560092112491126     (/mV2) 
  d1m = 4.0111125861771164e-06     (/mV3) 
  b2m = 0.02835498389099128     (/mV) 
  c2m = -0.00010981045750490155     (/mV2) 
  d2m = 7.972007758958973e-08     (/mV3) 
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