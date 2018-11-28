TITLE  squid sodium, potassium delayed rectifier, and potassium A channels
 
UNITS {
        (molar) = (1/liter)
        (S) = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
         F = (faraday) (coulomb)
         R = (mole k) (mV-coulomb/degC)
        (mM) =  (millimolar)
}
 
NEURON {
        SUFFIX hh3
        USEION na READ nao,nai WRITE ina
        USEION k READ ek WRITE ik
        RANGE  gnabar,gkhhbar,gkabar,ina,ikhh,ika,ik,ena,niv,miv,hiv,htv1,htv2 
        GLOBAL minf,hinf,ninf
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v   (mV)
        dt  (ms)
	:nai (mM)
	celsius = 35.0 (degC)
        gnabar  = 550.0e-6 (S/cm2)
        gkhhbar = 0.0 (S/cm2)
        gkabar  = 0.0  (S/cm2)
        :ek  = -90.0  (mV)
        :nao =  145  (mM)
        qv = 60.0 (mV)
        qs = 5.0  (1)
	  
 	
}
 
STATE {
        m <1e-4> h <1e-4> n <1e-4> p <1e-4> q <1e-4> 
}
 
ASSIGNED {
	ek (mV)
	nao (mM)
	nai (mM)
        ina (mA/cm2)
        ik (mA/cm2)
        ika (mA/cm2)
        ikhh (mA/cm2)
        ena (mV)
        minf hinf ninf qinf pinf
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        ena = R*(celsius+273.15)/F*log(nao/nai)
        ina = gnabar*m*m*m*h*(v - ena)
        ikhh = gkhhbar*n*n*n*(v - ek)      
        ika = gkabar*p*p*p*q*(v - ek)      
        ik = ika :+ ikhh
}
 
UNITSOFF
 
INITIAL {
        m = boltz(v,-34.6,6.0)
        h = boltz(v,-56.8,-7.8)
        n = boltz(v,-25.0,12.0)
        p = boltz(v,-43,24.0)
        q = boltz(v,-qv,-qs)
}

DERIVATIVE states {  :Computes state variables m, h, and n 
LOCAL minf,hinf,ninf,pinf,qinf,mtau,htau,ntau,ptau,qtau
        minf = boltz(v,-34.6,6.0)
        hinf = boltz(v,-56.8,-7.8)
        ninf = boltz(v,-25,12.0)
        pinf = boltz(v,-43,24.0)
        qinf = boltz(v,-qv,-qs)
        mtau = 1*boltz(v,-45.0,-1.5) - 1*boltz(v,-65.0,-0.5) +0.04
        htau = 3*56.0*boltz(v,-29.0,-4.5) - 3*56.0*boltz(v,-49.0,-2.0) +2.0
        ntau = 19.0/(1.0 + exp((v+ 39.0)/8.0)) - 0.0/(1.0 + exp((v+ 59.0)/20.0))+1.0
        ptau = 2.*exp(-(v+50)*(v+50)/550)+1.1
        qtau = 20.0
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
        p' = (pinf-p)/ptau
        q' = (qinf-q)/qtau
}
 
 
 
FUNCTION boltz(x,y,z) {
                boltz = 1/(1 + exp(-(x - y)/z))
}
 
UNITSON

