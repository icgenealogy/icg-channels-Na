NEURON
{
  SUFFIX HH 
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

  am = 0.20000115275005056     (/mV) 
  bm = -7.983777242388944     (1) 
  vhm = -43.44773455666261     (mV) 
  Am = 0.6051311698518747     (/ms) 
  b1m = -0.0796476210740026     (/mV) 
  c1m = 0.0007947737233582971     (/mV2) 
  d1m = -2.858554721413399e-06     (/mV3) 
  b2m = -0.10569489389695091     (/mV) 
  c2m = -0.0021318239184178536     (/mV2) 
  d2m = -1.7462681850189862e-05     (/mV3) 

  ah = -0.11331803737697743     (/mV) 
  bh = 7.634746070177136     (1) 
  vhh = -63.51426415607448     (mV) 
  Ah = 6.307860560314879     (/ms) 
  b1h = -0.05569848450001585     (/mV) 
  c1h = 0.00015502632040045064     (/mV2) 
  d1h = 3.1067882781038157e-07     (/mV3) 
  b2h = -0.051795271286430515     (/mV) 
  c2h = 0.0002737609742090115     (/mV2) 
  d2h = 9.298924082120821e-07     (/mV3) 

  an = 0.08857688150170923     (/mV) 
  bn = -2.907651190656495     (1) 
  vhn = -56.637661635469854     (mV) 
  An = 2.145156408698536     (/ms) 
  b1n = -0.03865928733255278     (/mV) 
  c1n = 0.00023514301366890335     (/mV2) 
  d1n = -6.053112261512405e-07     (/mV3) 
  b2n = -0.06934069968251609     (/mV) 
  c2n = -0.0012162047993050818     (/mV2) 
  d2n = -1.196498327880272e-05     (/mV3) 
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