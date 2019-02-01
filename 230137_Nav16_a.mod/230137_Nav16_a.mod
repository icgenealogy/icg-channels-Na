TITLE Nav1.6 ionic voltage-gated channel with kinetic scheme

COMMENT
A six-state markovian kinetic model of ionic channel.
Part of a study on kinetic models.
Author: Piero Balbi, August 2016
ENDCOMMENT

NEURON {
	SUFFIX na16a
	USEION na READ ena WRITE ina
	RANGE gbar, ina, g
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	ena (mV)
	celsius (degC)
	gbar  = 0.1	 (mho/cm2)
	
	C1C2b2	  = 14
	C1C2v2    = -8
	C1C2k2	  = -10
	
	C2C1b1	  = 2
	C2C1v1    = -38
	C2C1k1	  = 9
	C2C1b2	  = 14
	C2C1v2    = -8
	C2C1k2	  = -10

	C2O1b2	  = 14
	C2O1v2    = -18
	C2O1k2	  = -10
	
	O1C2b1	  = 4
	O1C2v1    = -48
	O1C2k1	  = 9
	O1C2b2	  = 14
	O1C2v2    = -18
	O1C2k2	  = -10
	
	C2O2b2	  = 0.0001
	C2O2v2	  = -10
	C2O2k2	  = -8
	
	O2C2b1	  = 0.0001
	O2C2v1	  = -55
	O2C2k1	  = 10
	O2C2b2	  = 0.0001
	O2C2v2	  = -20
	O2C2k2	  = -5
	
	O1I1b1	  = 6
	O1I1v1	  = -40
	O1I1k1	  = 13
	O1I1b2	  = 10
	O1I1v2	  = 15
	O1I1k2	  = -18
	
	I1O1b1	  = 0.00001
	I1O1v1	  = -40
	I1O1k1	  = 10
	
	I1C1b1	  = 0.1
	I1C1v1	  = -86
	I1C1k1	  = 9
	
	C1I1b2	  = 0.08
	C1I1v2	  = -55
	C1I1k2	  = -12
	
	I1I2b2	  = 0.00022
	I1I2v2	  = -50
	I1I2k2	  = -5

	I2I1b1	  = 0.0018
	I2I1v1	  = -90
	I2I1k1	  = 30
	
}

ASSIGNED {
	ina  (mA/cm2)
	g   (mho/cm2)
	
	C1C2_a (/ms)
	C2C1_a (/ms)
	C2O1_a (/ms)
	O1C2_a (/ms)
	C2O2_a (/ms)
	O2C2_a (/ms)
	O1I1_a (/ms)
	I1O1_a (/ms)
	I1I2_a (/ms)
	I2I1_a (/ms)
	I1C1_a (/ms)
	C1I1_a (/ms)
	
	Q10 (1)
}

STATE {
	C1
	C2
	O1
	O2
	I1
	I2
}


INITIAL {
	Q10 = 3^((celsius-20(degC))/10 (degC))
	SOLVE kin
	STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar * (O1 + O2)	: (mho/cm2)
	ina = g * (v - ena)   	: (mA/cm2)
}

KINETIC kin {
	rates(v)
	
	~ C1 <->  C2 (C1C2_a, C2C1_a)
	~ C2 <->  O1 (C2O1_a, O1C2_a)
	~ C2 <->  O2 (C2O2_a, O2C2_a)
	~ O1 <->  I1 (O1I1_a, I1O1_a)
	~ I1 <->  C1 (I1C1_a, C1I1_a)
	~ I1 <->  I2 (I1I2_a, I2I1_a)
	
	CONSERVE O1 + O2 + C1 + C2 + I1 + I2 = 1
}

FUNCTION rates2(v, b, vv, k) {
	rates2 = (b/(1+exp((v-vv)/k)))
}

PROCEDURE rates(v(mV)) {
UNITSOFF
	C1C2_a = Q10*(rates2(v, C1C2b2, C1C2v2, C1C2k2))
	C2C1_a = Q10*(rates2(v, C2C1b1, C2C1v1, C2C1k1) + rates2(v, C2C1b2, C2C1v2, C2C1k2))
	C2O1_a = Q10*(rates2(v, C2O1b2, C2O1v2, C2O1k2))
	O1C2_a = Q10*(rates2(v, O1C2b1, O1C2v1, O1C2k1) + rates2(v, O1C2b2, O1C2v2, O1C2k2))
	C2O2_a = Q10*(rates2(v, C2O2b2, C2O2v2, C2O2k2))
	O2C2_a = Q10*(rates2(v, O2C2b1, O2C2v1, O2C2k1) + rates2(v, O2C2b2, O2C2v2, O2C2k2))
	O1I1_a = Q10*(rates2(v, O1I1b1, O1I1v1, O1I1k1) + rates2(v, O1I1b2, O1I1v2, O1I1k2))
	I1O1_a = Q10*(rates2(v, I1O1b1, I1O1v1, I1O1k1))
	I1C1_a = Q10*(rates2(v, I1C1b1, I1C1v1, I1C1k1))
	C1I1_a = Q10*(rates2(v, C1I1b2, C1I1v2, C1I1k2))
	I1I2_a = Q10*(rates2(v, I1I2b2, I1I2v2, I1I2k2))
	I2I1_a = Q10*(rates2(v, I2I1b1, I2I1v1, I2I1k1))
UNITSON
}