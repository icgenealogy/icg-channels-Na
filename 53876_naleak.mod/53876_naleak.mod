COMMENT
This file, naleak.mod, implements the na leak current G_Na(leak) in 
Quadroni and Knopfel 1994 table 2
ENDCOMMENT

NEURON {
  SUFFIX naleak
  :NONSPECIFIC_CURRENT i
  USEION na READ ena WRITE ina
  RANGE g
}

PARAMETER {
  g = 30.2e-6 (siemens/cm2) < 0, 1e9 >
  :Erev = 50 (millivolt)
}

ASSIGNED {
  ena (millivolt)
  ina (milliamp/cm2)
  v (millivolt)
}

BREAKPOINT { ina = g * (v - ena) }
