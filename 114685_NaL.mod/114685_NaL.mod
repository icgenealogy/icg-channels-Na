TITLE sodium leak for GPi model neuron

COMMENT 
 
 Sodium leak.  This consists of a TTX-insensitive sodium current, which
 is active between spikes, has some voltage dependence, and underlies the
 resting membrane potential (AthertonBevan2005).  Other cation may be
 involved but for simplicity, only Na+ current was modeled.

ENDCOMMENT

UNITS {
   (mv) = (millivolt)
   (mA) = (milliamp)
}

NEURON {
    SUFFIX NaL
    USEION na READ ena,nai WRITE ina
    RANGE gna,inaL
    GLOBAL gbar_naL
}

PARAMETER {
    v (mV)
    gna = 2.8e-5   (mho/cm2)
    inaL = 0.0      (mA/cm2)
    ena
    nai
    celsius
}

ASSIGNED { 
    ina (mA/cm2)
    gbar_naL
}

BREAKPOINT {
    ina	= gna*gbar_naL*(v-ena)
    inaL = ina
}

INITIAL {
    gbar_naL = 1
}