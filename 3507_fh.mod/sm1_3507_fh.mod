NEURON
{
  SUFFIX fh 
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

  am = 0.10378543428739473     (/mV) 
  bm = -3.4864002446604005     (1) 
  vhm = -27.776669566278027     (mV) 
  Am = 0.15045721518934818     (/ms) 
  b1m = 0.9896486573318002     (/mV) 
  c1m = 0.037980766688067984     (/mV2) 
  d1m = 0.0003772291263004496     (/mV3) 
  b2m = 0.06780186943420956     (/mV) 
  c2m = -0.0009384859413805892     (/mV2) 
  d2m = 4.063954597941642e-06     (/mV3) 

  ap = 0.07498851579964118     (/mV) 
  bp = -0.9117188027125344     (1) 
  vhp = -22.258659353602596     (mV) 
  Ap = 1.1773289853272149     (/ms) 
  b1p = -0.02605747928519922     (/mV) 
  c1p = 0.00013623062773583279     (/mV2) 
  d1p = -3.350524115292254e-07     (/mV3) 
  b2p = -0.04734263945243467     (/mV) 
  c2p = -0.0001644790046953215     (/mV2) 
  d2p = 3.293517370642326e-08     (/mV3) 

  an = 0.11896137868177709     (/mV) 
  bn = -4.737712908380659     (1) 
  vhn = -42.19559986217699     (mV) 
  An = 0.9796366456155929     (/ms) 
  b1n = 0.07349133700118163     (/mV) 
  c1n = 0.0007133189615771017     (/mV2) 
  d1n = 3.3721331226036947e-06     (/mV3) 
  b2n = 0.04882128050208524     (/mV) 
  c2n = -0.00035373768258432827     (/mV2) 
  d2n = 1.0471958748721366e-06     (/mV3) 

  ah = -0.2148864772942064     (/mV) 
  bh = 13.48143328233753     (1) 
  vhh = -60.770584056182784     (mV) 
  Ah = 1.5374891297927484     (/ms) 
  b1h = -0.12369579620554436     (/mV) 
  c1h = 0.0012726012860188101     (/mV2) 
  d1h = -4.072733277901064e-06     (/mV3) 
  b2h = -0.10638945790944405     (/mV) 
  c2h = -0.0006336304039842476     (/mV2) 
  d2h = 3.979873533435403e-06     (/mV3) 
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
  pInf 
  pTau 
  nInf 
  nTau 
  hInf 
  hTau 
}

STATE
{
  m
  p
  n
  h
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m*p*p*n*n*h
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
  p' = (pInf - p) / pTau 
  n' = (nInf - n) / nTau 
  h' = (hInf - h) / hTau 
}

INITIAL
{
  rates(v)
  m = mInf 
  p = pInf 
  n = nInf 
  h = hInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    pInf = 1/(1 + exp(-ap*v + bp)) 
    pTau = Ap / ( exp(-(b1p*(v-vhp) + c1p*(v-vhp)^2 + d1p*(v-vhp)^3)) + exp((b2p*(v-vhp) + c2p*(v-vhp)^2 + d2p*(v-vhp)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 


  UNITSON
}