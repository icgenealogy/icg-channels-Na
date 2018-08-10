NEURON
{
  SUFFIX GrG_Nar 
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

  as = 0.05335578126819678     (/mV) 
  bs = 0.7669185828500686     (1) 
  vhs = -49.07286492411861     (mV) 
  As = 2.8926794627241743     (/ms) 
  b1s = 0.017442135653854378     (/mV) 
  c1s = -0.0006190383102108388     (/mV2) 
  d1s = 2.755537512346411e-06     (/mV3) 
  b2s = -0.26186989704339847     (/mV) 
  c2s = -0.011022637319068123     (/mV2) 
  d2s = -0.00013942860931915493     (/mV3) 

  af = -0.07828373239906448     (/mV) 
  bf = 3.0236323280518014     (1) 
  vhf = -25.386294924211782     (mV) 
  Af = 0.6550689141916181     (/ms) 
  b1f = 0.00295002486272643     (/mV) 
  c1f = -4.811593364764066e-05     (/mV2) 
  d1f = -2.2879081598262235e-07     (/mV3) 
  b2f = 0.11766252775559673     (/mV) 
  c2f = 1.6645268149531263e-05     (/mV2) 
  d2f = -1.9092427998244067e-05     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  ena	(mV)
  ina	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  sInf 
  sTau 
  fInf 
  fTau 
}

STATE
{
  s
  f
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*s*f
  ina = g*(v-ena)
}

DERIVATIVE states
{
  rates(v)
  s' = (sInf - s) / sTau 
  f' = (fInf - f) / fTau 
}

INITIAL
{
  rates(v)
  s = sInf 
  f = fInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    sInf = 1/(1 + exp(-as*v + bs)) 
    sTau = As / ( exp(-(b1s*(v-vhs) + c1s*(v-vhs)^2 + d1s*(v-vhs)^3)) + exp((b2s*(v-vhs) + c2s*(v-vhs)^2 + d2s*(v-vhs)^3)) ) 

    fInf = 1/(1 + exp(-af*v + bf)) 
    fTau = Af / ( exp(-(b1f*(v-vhf) + c1f*(v-vhf)^2 + d1f*(v-vhf)^3)) + exp((b2f*(v-vhf) + c2f*(v-vhf)^2 + d2f*(v-vhf)^3)) ) 


  UNITSON
}