NEURON {
    SUFFIX na }
NEURON {
    USEION na READ ena WRITE ina }
ASSIGNED {
    ena (mV)
    ina }
PARAMETER {
        cai (mM)
        cao (mM)
        v (mV)
        celsius  = 37 (degC)
        dt (ms)
	:erev 		= 45        (mV)
	gmax 		= 0.1     (umho)
        vrest           = 0

	mvalence 	= 4.3
	mgamma 		=  0.7
	mbaserate 	=  4.2
	mvhalf 		=  -38
	mbasetau 	=  0.05
	mtemp 		=  37
        mq10            =  3
	mexp 		=  3

	hvalence 	= -6
	hgamma		=  0.5
	hbaserate 	=  0.2
	hvhalf 		=  -42
	hbasetau 	=  0.5
	htemp 		=  37
        hq10            =  3
	hexp 		=  1



	vmax 		= 50     (mV)
	vmin 		= -100   (mV)
} : end PARAMETER


INCLUDE "custom_code/inc_files/64229_bg_cvode.inc"
PROCEDURE iassign() { i = g * (v - ena) ina=i }

:** nap
