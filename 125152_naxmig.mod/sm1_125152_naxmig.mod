NEURON
{
  SUFFIX naxMig 
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

  ah = -0.2499995533103707     (/mV) 
  bh = 12.4999756735426     (1) 
  vhh = -48.588211938323816     (mV) 
  Ah = 15.020379746520636     (/ms) 
  b1h = -0.19841291761395918     (/mV) 
  c1h = 0.002760874137099447     (/mV2) 
  d1h = -1.103018246587544e-05     (/mV3) 
  b2h = -0.22579872724729372     (/mV) 
  c2h = -0.0071527961627621435     (/mV2) 
  d2h = -7.864845094713633e-05     (/mV3) 

  am = 0.138888734047016     (/mV) 
  bm = -4.2267264564045925     (1) 
  vhm = -29.969687497478446     (mV) 
  Am = 0.2734937739472201     (/ms) 
  b1m = 0.05003030632720269     (/mV) 
  c1m = 0.0004959323747667389     (/mV2) 
  d1m = 1.9052334755346798e-06     (/mV3) 
  b2m = 0.08065772501695206     (/mV) 
  c2m = -0.0008879448561441898     (/mV2) 
  d2m = 3.089437660196013e-06     (/mV3) 
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