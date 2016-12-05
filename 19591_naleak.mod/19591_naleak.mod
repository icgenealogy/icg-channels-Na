TITLE Ca leak

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT { v FROM -100 TO 50 WITH 50	(mV) }

NEURON {
	SUFFIX naleak
	:NONSPECIFIC_CURRENT i
	USEION na READ ena WRITE ina
        RANGE gbar
}

PARAMETER {
	gbar = 1.167e-5	(mho/cm2)
	:e = 55	(mV)
}

ASSIGNED {
    ina    (mA/cm2)
    ena (mV)
}

BREAKPOINT {
	ina = gbar*(v - ena)
}
