TITLE inad :Transient sodium current in the interneuron dendrites

NEURON {
	SUFFIX inad
	USEION na READ ena WRITE ina
	RANGE gna, ina, qna
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
        PI   = (pi) (1)
	FARADAY	= 96485.309 (coul/mole)
}

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {
	celsius=37	(degC)
	gna=1e-3	(mho/cm2)
}

ASSIGNED { 
	ina	(mA/cm2)
	v	(mV)
	ena	(mV)
        diam	(um)
}

STATE { ma mb ha hb qna }

BREAKPOINT {
	SOLVE nastate METHOD sparse
	ina = gna*ma*ma*ha*(v-ena)
}

INITIAL {
	ma=m_inf(v)
	ha=h_inf(v)
	mb=1-ma
	hb=1-ha
        qna = 0
	ina = gna*ma*ma*ha*(v-ena)
}

LOCAL a1, a2, b1, b2

KINETIC nastate {	
	a1 = m_a(v)
	a2 = m_b(v)
	b1 = h_a(v)
	b2 = h_b(v)
	~ mb <-> ma (a1, a2)
	~ hb <-> ha (b1, b2)
	CONSERVE ma + mb = 1
	CONSERVE ha + hb = 1

	COMPARTMENT diam*diam*PI/4 { qna }
        ~ qna << (-ina*PI*diam*(1e4)/FARADAY)
}

FUNCTION m_a(v(mV)) {
	TABLE FROM -150 TO 150 WITH 301
	m_a=0.32*(13.1-v-60)/(exp((13.1-v-60)/4)-1)
}

FUNCTION m_b(v(mV)) {
	TABLE FROM -150 TO 150 WITH 301
	m_b=0.28*(v-40.1+60)/(exp((v-40.1+60)/5)-1)
}

FUNCTION h_a(v(mV)) {
	TABLE FROM -150 TO 150 WITH 301
	h_a = 0.128*exp((17-v-60)/18)
}

FUNCTION h_b(v(mV)) {
	TABLE FROM -150 TO 150 WITH 301
	h_b = 4/(1+exp((40-v-60)/5))
}

FUNCTION m_inf(v(mV)) {
	m_inf = m_a(v)/(m_a(v)+m_b(v))
}

FUNCTION h_inf(v(mV)) {
	h_inf = h_a(v)/(h_a(v)+h_b(v))
}
