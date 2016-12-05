TITLE Fast Sodium current 
: Unit check passed

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 
NEURON { 
	SUFFIX NaF
	USEION na READ ena WRITE ina
	RANGE g, gmax, ina
	GLOBAL minf, mtau, hinf, htau
}
PARAMETER { : PARAMETERS ARE BY DEFAULT GLOBAL VARIABLES
	gmax = 0.0195 	(mho/cm2)
	ena (mV)
	m_vh = -23.9	(mV)	: half activation 
	m_ve = -11.8	(mV)	: slope
	h_vh = -62.9	(mV) : half activation
	h_ve = 10.7	(mV) : slope
} 
ASSIGNED { 
	v 		(mV)
	g		(mho/cm2)
	ina 		(mA/cm2) 
	minf 		(1)
	mtau 	(ms)
	hinf 		(1)
	htau 		(ms) 
} 
STATE {
	m h
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	g = gmax * m*m*m * h
	ina = g * (v - ena ) 
} 
INITIAL { 
	rates(v) 
	m = minf
	h = hinf
} 
DERIVATIVE states { 
	rates(v)
	m' = ( minf - m ) / mtau 
	h' = ( hinf - h ) / htau
}
FUNCTION_TABLE tabmtau(v(mV)) (ms) 
FUNCTION_TABLE tabhtau(v(mV)) (ms) 

PROCEDURE rates(v(mV)) { 
	TABLE mtau, htau, minf, hinf DEPEND h_vh FROM -120 TO 40 WITH 160
:	TABLE mtau, htau, minf, hinf DEPEND h_vh FROM -120 TO 40 WITH 320
	mtau = tabmtau(v)
	htau = tabhtau(v)
	minf = 1/(1 + exp((v - m_vh)/m_ve))
	hinf = 1/(1 + exp((v - h_vh)/h_ve))
}
