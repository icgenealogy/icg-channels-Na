COMMENT

   **************************************************
   File generated by: neuroConstruct v1.7.1 
   **************************************************

   This file holds the implementation in NEURON of the Cell Mechanism:
   NaP_iAMC_Fig10Hii_ChannelML (Type: Channel mechanism, Model: ChannelML based process)

   with parameters: 
   /channelml/@units = Physiological Units 
   /channelml/notes = Mitral Cell Persistent Sodium ion Channel 
   /channelml/channel_type/@name = NaP_iAMC_Fig10Hii_ChannelML 
   /channelml/channel_type/@density = yes 
   /channelml/channel_type/status/@value = stable 
   /channelml/channel_type/status/comment = Sodium Persistent conductance modified from Rubin an Cleland 2006 using AOB mitral cell data from the Marc Spehr RWTH Aachen 
   /channelml/channel_type/status/contributor/name = Simon O'Connor 
   /channelml/channel_type/notes = Na Channel 
   /channelml/channel_type/authorList/modelTranslator/name = Simon O'Connor 
   /channelml/channel_type/authorList/modelTranslator/institution = UH 
   /channelml/channel_type/authorList/modelTranslator/email = simon.oconnor - at - btinternet.com 
   /channelml/channel_type/current_voltage_relation/@cond_law = ohmic 
   /channelml/channel_type/current_voltage_relation/@ion = na 
   /channelml/channel_type/current_voltage_relation/@default_gmax = 0.06 
   /channelml/channel_type/current_voltage_relation/@default_erev = 67 
   /channelml/channel_type/current_voltage_relation/@charge = 1 
   /channelml/channel_type/current_voltage_relation/gate[1]/@name = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/@instances = 3 
   /channelml/channel_type/current_voltage_relation/gate[1]/closed_state/@id = m0 
   /channelml/channel_type/current_voltage_relation/gate[1]/open_state/@id = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/open_state/@fraction = 1 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@name = tau 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@from = m0 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@to = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@expr_form = generic 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@expr = (1+(4 * (exp(0 - ((v + 50)/20)^2)))) 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@name = inf 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@from = m0 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@to = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@expr_form = sigmoid 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@rate = 0.499622025796 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@scale = -4.9 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@midpoint = -59.0 
   /channelml/channel_type/current_voltage_relation/gate[2]/@name = h 
   /channelml/channel_type/current_voltage_relation/gate[2]/@instances = 1 
   /channelml/channel_type/current_voltage_relation/gate[2]/closed_state/@id = h0 
   /channelml/channel_type/current_voltage_relation/gate[2]/open_state/@id = h 
   /channelml/channel_type/current_voltage_relation/gate[2]/open_state/@fraction = 1 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@name = tau 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@from = h0 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@to = h 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@expr_form = generic 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@expr = (5000+(16000 * (exp(0 - ((v + 50)/20)^2)))) 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@name = inf 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@from = h0 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@to = h 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@expr_form = sigmoid 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@rate = 0.499622025796 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@scale = 4.9 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@midpoint = -59.0 
   /channelml/channel_type/current_voltage_relation/gate[3]/@name = n 
   /channelml/channel_type/current_voltage_relation/gate[3]/@instances = 1 
   /channelml/channel_type/current_voltage_relation/gate[3]/closed_state/@id = n0 
   /channelml/channel_type/current_voltage_relation/gate[3]/open_state/@id = n 
   /channelml/channel_type/current_voltage_relation/gate[3]/open_state/@fraction = 1 
   /channelml/channel_type/current_voltage_relation/gate[3]/time_course/@name = tau 
   /channelml/channel_type/current_voltage_relation/gate[3]/time_course/@from = n0 
   /channelml/channel_type/current_voltage_relation/gate[3]/time_course/@to = n 
   /channelml/channel_type/current_voltage_relation/gate[3]/time_course/@expr_form = generic 
   /channelml/channel_type/current_voltage_relation/gate[3]/time_course/@expr = (2+(4 * (exp(0 - ((v + 50)/20)^2)))) 
   /channelml/channel_type/current_voltage_relation/gate[3]/steady_state/@name = inf 
   /channelml/channel_type/current_voltage_relation/gate[3]/steady_state/@from = n0 
   /channelml/channel_type/current_voltage_relation/gate[3]/steady_state/@to = n 
   /channelml/channel_type/current_voltage_relation/gate[3]/steady_state/@expr_form = sigmoid 
   /channelml/channel_type/current_voltage_relation/gate[3]/steady_state/@rate = 0.499622025796 
   /channelml/channel_type/current_voltage_relation/gate[3]/steady_state/@scale = 4.9 
   /channelml/channel_type/current_voltage_relation/gate[3]/steady_state/@midpoint = -59.0 
   /channelml/channel_type/impl_prefs/table_settings/@max_v = 100 
   /channelml/channel_type/impl_prefs/table_settings/@min_v = -100 
   /channelml/channel_type/impl_prefs/table_settings/@table_divisions = 2000 

// File from which this was generated: /home/Simon/NML2_Test/iAMC_Fig10H2T/AOB_MC_neuroConstruct/cellMechanisms/NaP_iAMC_Fig10Hii_ChannelML/NaChannel.xml

// XSL file with mapping to simulator: /home/Simon/NML2_Test/iAMC_Fig10H2T/AOB_MC_neuroConstruct/cellMechanisms/NaP_iAMC_Fig10Hii_ChannelML/ChannelML_v1.8.1_NEURONmod.xsl

ENDCOMMENT


?  This is a NEURON mod file generated from a ChannelML file

?  Unit system of original ChannelML file: Physiological Units

COMMENT
    Mitral Cell Persistent Sodium ion Channel
ENDCOMMENT

TITLE Channel: NaP_iAMC_Fig10Hii_ChannelML

COMMENT
    Na Channel
ENDCOMMENT


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (um) = (micrometer)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (l) = (liter)
}


    
NEURON {
      

    SUFFIX NaP_iAMC_Fig10Hii_ChannelML
    USEION na READ ena WRITE ina VALENCE 1 ? reversal potential of ion is read, outgoing current is written
           
        
    RANGE gmax, gion
    
    RANGE minf, mtau
    
    RANGE hinf, htau
    
    RANGE ninf, ntau
    
}

PARAMETER { 
      

    gmax = 0.000059999999999999995 (S/cm2)  ? default value, should be overwritten when conductance placed on cell
    
}



ASSIGNED {
      

    v (mV)
    
    celsius (degC)
          

    ? Reversal potential of na
    ena (mV)
    ? The outward flow of ion: na calculated by rate equations...
    ina (mA/cm2)
    
    
    gion (S/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
    ninf
    ntau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
        
    gion = gmax * ((1*m)
^3) * ((1*h)
^1) * ((1*n)
^1)      

    ina = gion*(v - ena)
            

}



INITIAL {
    
    ena = 67
        
    rates(v)
    m = minf
        h = hinf
        n = ninf
        
    
}
    
STATE {
    m
    h
    n
    
}



DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
            h' = (hinf - h)/htau
            n' = (ninf - n)/ntau
            

}

PROCEDURE rates(v(mV)) {  
    
    ? Note: not all of these may be used, depending on the form of rate equations
    LOCAL  alpha, beta, tau, inf, gamma, zeta
, temp_adj_m,
         A_inf_m, B_inf_m, Vhalf_inf_m
, temp_adj_h,
         A_inf_h, B_inf_h, Vhalf_inf_h
, temp_adj_n,
         A_inf_n, B_inf_n, Vhalf_inf_n
    
    TABLE minf, mtau,hinf, htau,ninf, ntau
 DEPEND celsius FROM -100 TO 100 WITH 2000
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_h = 1
    temp_adj_n = 1
    
            
                
           

        
    ?      ***  Adding rate equations for gate: m  ***
         
    ? Found a generic form of the rate equation for tau, using expression: (1+(4 * (exp(0 - ((v + 50)/20)^2))))
    tau = (1+(4 * (exp(0 - ((v + 50)/20)^2))))
        
    mtau = tau/temp_adj_m
    
    ? Found a parameterised form of rate equation for inf, using expression: A / (1 + exp((v-Vhalf)/B))
    A_inf_m = 0.499622025796
    B_inf_m = -4.9
    Vhalf_inf_m = -59.0 
    inf = A_inf_m / (exp((v - Vhalf_inf_m) / B_inf_m) + 1)
    
    minf = inf
    


    ?     *** Finished rate equations for gate: m ***
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate: h  ***
         
    ? Found a generic form of the rate equation for tau, using expression: (5000+(16000 * (exp(0 - ((v + 50)/20)^2))))
    tau = (5000+(16000 * (exp(0 - ((v + 50)/20)^2))))
        
    htau = tau/temp_adj_h
    
    ? Found a parameterised form of rate equation for inf, using expression: A / (1 + exp((v-Vhalf)/B))
    A_inf_h = 0.499622025796
    B_inf_h = 4.9
    Vhalf_inf_h = -59.0 
    inf = A_inf_h / (exp((v - Vhalf_inf_h) / B_inf_h) + 1)
    
    hinf = inf
    


    ?     *** Finished rate equations for gate: h ***
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate: n  ***
         
    ? Found a generic form of the rate equation for tau, using expression: (2+(4 * (exp(0 - ((v + 50)/20)^2))))
    tau = (2+(4 * (exp(0 - ((v + 50)/20)^2))))
        
    ntau = tau/temp_adj_n
    
    ? Found a parameterised form of rate equation for inf, using expression: A / (1 + exp((v-Vhalf)/B))
    A_inf_n = 0.499622025796
    B_inf_n = 4.9
    Vhalf_inf_n = -59.0 
    inf = A_inf_n / (exp((v - Vhalf_inf_n) / B_inf_n) + 1)
    
    ninf = inf
    


    ?     *** Finished rate equations for gate: n ***
    

         

}


UNITSON


