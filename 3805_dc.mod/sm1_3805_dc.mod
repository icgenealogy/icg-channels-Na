NEURON
{
  SUFFIX dc 
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

  ah = -0.14147646559614385     (/mV) 
  bh = -1.8166399762589642     (1) 
  vhh = 16.970913633245225     (mV) 
  Ah = 3.5637056559164106     (/ms) 
  b1h = -0.11489572940349281     (/mV) 
  c1h = 0.001623578610017657     (/mV2) 
  d1h = -6.930964869364258e-06     (/mV3) 
  b2h = -0.029527254594561864     (/mV) 
  c2h = 0.0003731628683317974     (/mV2) 
  d2h = 2.1832741203014354e-06     (/mV3) 

  an = 0.057782437320304025     (/mV) 
  bn = 1.3851474355495983     (1) 
  vhn = 25.843991211390314     (mV) 
  An = 4.314217290092606     (/ms) 
  b1n = 0.00024977933288029597     (/mV) 
  c1n = -0.00022017292054386176     (/mV2) 
  d1n = -1.1587133830381694e-06     (/mV3) 
  b2n = 0.03862975254674466     (/mV) 
  c2n = -0.0001253717248449498     (/mV2) 
  d2n = -3.413417681727764e-06     (/mV3) 

  am = 0.1430880355709965     (/mV) 
  bm = 5.008094809338396     (1) 
  vhm = 41.12745384959255     (mV) 
  Am = 0.21370108358860276     (/ms) 
  b1m = -0.06903637750088314     (/mV) 
  c1m = 0.001118555754038185     (/mV2) 
  d1m = -7.509566502774477e-06     (/mV3) 
  b2m = -0.04532455444121375     (/mV) 
  c2m = -0.0003947201577052832     (/mV2) 
  d2m = -1.373020114180471e-06     (/mV3) 
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
  nInf 
  nTau 
  mInf 
  mTau 
}

STATE
{
  h
  n
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*n*n*n*n*m*m*m
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  n' = (nInf - n) / nTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  n = nInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}