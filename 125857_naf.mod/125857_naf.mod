COMMENT

	Maciej T. Lazarewicz, mlazarew@seas.upenn.edu

ENDCOMMENT

NEURON {

	SUFFIX naf
	USEION na READ ena WRITE ina
	RANGE gnabar, ina
}
	
UNITS {

	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {

    gnabar  = 30 (mS/cm2)
    :ena  = 55 (mV)
}
    
ASSIGNED {
    ena (mV)
    v    (mV)
    ina  (mA/cm2)
    minf (1)
    hinf (1)
    tauh (ms)
}

STATE { h }

INITIAL {
    
    rates(v)
    h  = hinf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ina = (1e-3) * gnabar * minf^2 * h * (v-ena)
}


DERIVATIVE states { 

    rates(v)
    h' = (hinf-h)/tauh
}


:ina
PROCEDURE rates(v(mV)) { LOCAL a, b

    a    = fun3(v,  -46.9, -0.32,    -4) 
    b    = fun3(v,  -19.9,  0.28,     5) 
    minf = a/(a+b)
    
    a    = fun1(v,  -43,    0.128,  -18) 
    b    = fun2(v,  -20,    4,       -5)
    hinf = a/(a+b)
    tauh = 1.0/(a+b)
}

COMMENT

	Maciej T. Lazarewicz, mlazarew@seas.upenn.edu

ENDCOMMENT



:-------------------------------------------------------------------
FUNCTION fun1(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

	 fun1 = A*exp((v-V0)/B)
}

FUNCTION fun2(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

	 fun2 = A/(exp((v-V0)/B)+1)
}

FUNCTION fun3(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

    if(fabs((v-V0)/B)<1e-6) {
    :if(v==V0) {
        fun3 = A*B/1(mV) * (1- 0.5 * (v-V0)/B)
    } else {
        fun3 = A/1(mV)*(v-V0)/(exp((v-V0)/B)-1)
    }
}

FUNCTION min(x,y) { if (x<=y){ min = x }else{ min = y } }
