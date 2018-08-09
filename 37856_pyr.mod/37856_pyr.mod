TITLE pyr.mod   DCN pyramidal cell model  
 
COMMENT

Revised version of DCN Pyramidal cell model based on new hh.hoc file in NEURON

This model implements a Dorsal Cochlear Nucleus Pyramidal point cell
based on kinetic data from Kanold and Manis (1999) and Kanold's dissertation (1999)

-- 15 Jan 1999 P. Manis

Added export of start states for some variables to do perturbation tests
These start values replace the "inf" values used in the initialization procedure
Note that if the start variable is set to a value less than 0, 
then the default initialization will be done. Typically I use a value of -1 for this flagging
Note also that it is possible to set the initial values > 1 but this is meaningless in terms of
the present equations. 
-- 5 Feb 1999 P. Manis

Removed slow Ih current 30 Jan 2000. P. Manis
- also renamed variables to saner forms
Added Patrick's version of ih as ihd

Added scale factor for rate constants. rsf = 1 means default values; rf > 1 slows things down
5/30/2000 P. Manis

ENDCOMMENT
 
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
? interface
NEURON {
        SUFFIX pyr
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
		USEION h READ eh WRITE ih VALENCE 1
		NONSPECIFIC_CURRENT il

		RANGE gna, gk, minf, hinf, ninf		: sodium channels and delayed rectifier
		GLOBAL mtau, htau, ntau, gbar, gkbar	: time constants for sodium channels and delayed rectifier

        RANGE gkis, kis_a_inf, kis_i_inf	: slow inactivating potassium current
		RANGE gkif, kif_a_inf, kif_i_inf	: fast inactivating potassium current
		RANGE akif, akis					: state variables, fast and slow
		GLOBAL kif_a_tau, kif_i_tau, gkifbar	: fast inactivating potassium current
		GLOBAL kis_a_tau, kis_i_tau, gkisbar	: slow inactivating potassium current
		GLOBAL kif_a_start, kif_i_start
		GLOBAL kis_a_start, kis_i_start 		:flags for perturbation analysis
		GLOBAL kif_vh
		
		RANGE gh, kh_m_inf, kh_n_inf, aih
		GLOBAL kh_m_tau, kh_n_tau, ghbar

		GLOBAL gl, el							: leak (linear, ohmic with reversal at el)
		
		GLOBAL rsf			: rate scale factor for fast transient current
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius (degC)
        dt (ms)
		ek (mV)
        ena (mV)
        gbar = 0.02857 (mho/cm2)	<0,1e9>
        gkbar = 0.0 (mho/cm2)	<0,1e9>
		gkifbar = 0.0 (mho/cm2) <0,1e9>
        gkisbar = 0.0 (mho/cm2)
        ghbar = 0.0 (mho/cm2) <0,1e9>
		eh (mV)
		gl = 0.0 (mho/cm2)	<0,1e9>
        el = -70 (mV) 
		mtau = 0.05 (ms) <0.01,100>
		htau = 0.5 (ms) <0.1,100>
		ntau = 0.5 (ms) <0.1,100>
		kif_vh = -89.6
		kif_a_start = -1
		kif_i_start = -1
		kis_a_start = -1
		kis_i_start = -1
	    rsf = 1
}
 
STATE {
        m h n kifa kifi kisa kisi khm khn
}
 
ASSIGNED {
	gna (mho/cm2)
	gk (mho/cm2)
	gkif (mho/cm2)
	gkis (mho/cm2)
	gh (mho/cm2)
    ina (mA/cm2)
    ik (mA/cm2)
	ih (mA/cm2)
    il (mA/cm2)
    minf hinf
	ninf 
	kif_a_inf kif_i_inf
	kis_a_inf kis_i_inf
	kif_a_tau kif_i_tau
	kis_a_tau kis_i_tau
	kh_m_inf kh_n_inf
	kh_m_tau kh_n_tau	
	akif
	akis
	aih
}

LOCAL mexp, hexp, nexp, kif_a_exp, kif_i_exp, kis_a_exp, kis_i_exp, kh_m_exp, kh_n_exp
 
? currents
BREAKPOINT {
    SOLVE states METHOD cnexp
    	gna = gbar*m*m*h
	ina = gna*(v - ena)
   	gk = gkbar*n*n
	akif = kifa*kifa*kifa*kifa*kifi
	gkif = gkifbar*akif
	akis = kisa*kisa*kisa*kisa*kisi
	gkis = gkisbar*akis
	aih = khm*khn
	gh = ghbar*aih
	ik = gk*(v - ek) + gkif*(v - ek) + gkis*(v - ek)
 	ih = gh*(v - eh)
	
	il = gl*(v - el)
}
? currents

UNITSOFF 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	khm = kh_m_inf
	khn = kh_n_inf
	if(kif_a_start < 0) {		: if xx_(i/a)_start is > 0, then perturbation is done at onset of computations.
		kifa = kif_a_inf
	} else {
		kifa = kif_a_start
	}
	if(kif_i_start < 0) {
		kifi = kif_i_inf
	} else {
		kifi = kif_i_start
	}
	if(kis_a_start < 0) {
		kisa = kis_a_inf
	} else {
		kisa = kis_a_start
	}
	if(kis_i_start < 0) {
		kisi = kis_i_inf
	} else {
		kisi = kis_i_start
	}
}


? states
DERIVATIVE states {  
    rates(v)
    m' = (minf - m) / mtau
    h' = (hinf - h) / htau
    n' = (ninf - n) / ntau
	kifa' = (kif_a_inf - kifa) / kif_a_tau*rsf
	kifi' = (kif_i_inf - kifi) / kif_i_tau*rsf
	kisa' = (kis_a_inf - kisa) / kis_a_tau
	kisi' = (kis_i_inf - kisi) / kis_i_tau
	khm' = (kh_m_inf - khm) / kh_m_tau
	khn' = (kh_n_inf - khn) / kh_n_tau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
	LOCAL  alpha, beta, sum
	TABLE minf, mtau, hinf, htau, ninf, ntau, kif_a_inf, kif_a_tau, kif_i_inf, kif_i_tau, kis_a_inf, kis_a_tau, kis_i_inf, kis_i_tau, kh_m_inf, kh_n_inf, kh_m_tau, kh_n_tau DEPEND celsius, kif_vh FROM -200 TO 100 WITH 400

UNITSOFF
	q10 = 3^((celsius - 32)/10)

	:"m" sodium activation system
		minf = na_m(v)


	:"h" sodium inactivation system
		hinf = na_h(v)

	:"n" potassium activation system
        ninf = kd_m(v)

	:"kif" fast inactivation potassium channel - activation and inactivation 
		kif_a_inf = kif_m(v)
		kif_i_inf = kif_h(v)
		kif_a_tau = kif_mt(v)
		kif_i_tau = kif_ht(v) 

	:"kis" slow inactivating potassium channel - activation and inactivation
		kis_a_inf = kis_m(v)
		kis_i_inf = kis_h(v)
		kis_a_tau = kis_mt(v)
		kis_i_tau = kis_ht(v)

	:"kh" adaptation of Destexhe hyp-activated cation current by Patrick Kanold
		kh_m_inf = kh_m(v) 
		kh_n_inf = kh_n(v)
		kh_m_tau = kh_mt(v)
		kh_n_tau = kh_nt(v) 
}
 
: Make these as functions so we can view them from hoc, although this 
: may slow things down a bit

FUNCTION na_m(x) { : sodium activation
	na_m = 1/(1+exp(-(x+38)/3.0))	: flat time constants
}

FUNCTION na_h(x) { : sodium inactivation
	na_h = 1/(1+exp((x+43)/3.0))	: flat time constants
}

FUNCTION kd_m(x) { : potassium inactivation
	kd_m = 1/(1+exp(-(x+40)/3))		: flat time constants
}

FUNCTION kif_m(x) { : ikif activation
	kif_m = 1/(1+exp(-(x+53)/25.8))
}	

FUNCTION kif_h(x) { : ikif inactivation
	kif_h = 1/(1+exp((x-kif_vh)/6.7))
}	

FUNCTION kif_mt(x) { : ikif activation tau
	kif_mt = 0.15*exp((x+57)/10) + 0.3*exp(-(x+57)/10)
	kif_mt = 0.5 + (1/kif_mt)
}

FUNCTION kif_ht(x) { : ikif inactivation tau
	kif_ht = 0.015*exp((x+87)/20)+0.03*exp(-(x+87)/20)
	kif_ht = 10 + (1/kif_ht)
}

FUNCTION kis_m(x) { : ikis activation
	kis_m = 1/(1+exp(-(x+40.9)/23.7))
}

FUNCTION kis_h(x) { : ikis inactivation
	kis_h = 1/(1+exp((x+38.4)/9))
}

FUNCTION kis_mt(x) { : ikis activation tau
	kis_mt = 0.15*exp((x+40)/10) + 0.3*exp(-(x+40)/10)
	kis_mt = 0.5 + 1/kis_mt
}

FUNCTION kis_ht(x) { : ikis inactivation tau
	kis_ht = 200
}

FUNCTION kh_m(x) {
	kh_m = 1/(1+exp((x+68.9)/6.5))
}

FUNCTION kh_n(x) {
	kh_n = 1/(1+exp((x+68.9)/6.5)) : same as kh_m, but for completeness, compute this
}

FUNCTION kh_mt(v) {
	kh_mt = exp((v+183.6)/15.24)
}

FUNCTION kh_nt(v) {
	kh_nt = exp((v+158.6)/11.2)/(1+exp((v+75)/5.5))
}


FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
