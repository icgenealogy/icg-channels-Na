UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mS) = (millisiemens)
}

NEURON {
  SUFFIX Nafbwb
  USEION na READ ena WRITE ina
  RANGE phih
  RANGE gna, ena, taoh : testing
  GLOBAL ena
}
	
PARAMETER {
  gna  = 35 (mS/cm2)
  ena  = 55 (mV)
  phih = 5
}
    
ASSIGNED {
  v       (mV)
  ina     (mA/cm2)
  minf    (1)
  hinf    (1)
  taoh    (ms)
  celsius (degC)
}

STATE { h }

PROCEDURE iassign () { ina = (1e-3) * gna * minf^3 * h * (v-ena) }

INITIAL {
  rates(v)
  h = hinf
  iassign()
}

BREAKPOINT {
  SOLVE states METHOD cnexp	
  iassign()
}

DERIVATIVE states { 
  rates(v)
  h' = (hinf-h)/taoh
}

PROCEDURE rates(v(mV)) { LOCAL am, bm, ah, bh, q10
    
  q10  = phih:^((celsius-27.0(degC))/10.0(degC))	
    
  am   = fun3(v,  -35, -0.1,    -10)
  bm   = fun1(v,  -60,  4,      -18) 
  minf = am/(am+bm)
 
  ah   = fun1(v,  -58,    0.07,  -20)
  bh   = fun2(v,  -28,    1,     -10)
  hinf = ah/(ah+bh)
  taoh = 1./((ah+bh)*q10)
}

INCLUDE "custom_code/inc_files/185858_aux_fun.inc"
