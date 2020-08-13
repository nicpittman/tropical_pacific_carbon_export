#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 10:04:45 2020
@author: Nic.Pittman@utas.edu.au

Python implementation of Wanninkhof (2014) CO2 flux measurements
Calculates Schmidt number, CO2 solubility and thus the approximate air-sea CO2 flux. 

Check the papers to convert for different gas coefficients. 

Be careful of units during the final F=k*Ko*DpCO2
    
    k1=(k*24*365)/100                    #cm/hr to m/yr              #Conversion = *24*365/100
    ko1=ko*1000                          #mol/atm/l to mol/atm/m3    #Conversion = *1000
    dpco21=dpco2 /1000000                #uatm to atm                #Conversion = *1 000 000
    
Example:
    
    print(carbon_flux(33,20,6,425,400))
where:
    carbon_flux(salinity,temperature (C),windspeed U,pco2atm,pco2sw,dpco2 (optional))
    
Source:
Sutton, A. J., Wanninkhof, R., Sabine, C. L., Feely, R. A., Cronin, M. F., & Weller, R. A. (2017). Variability and trends in surface seawater p CO 2 and CO 2 flux in the Pacific Ocean: Pacific Ocean CO 2 Flux Trends. Geophysical Research Letters, 44(11), 5627–5636. https://doi.org/10.1002/2017GL073814
Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean revisited. Limnology and Oceanography: Methods, 12(6), 351–362. https://doi.org/10.4319/lom.2014.12.351
Weiss, R. F. (1974). Carbon dioxide in water and seawater: The solubility of a non-ideal gas. Marine Chemistry, 2(3), 203–215. https://doi.org/10.1016/0304-4203(74)90015-2
"""
import numpy as np

def schmidt(t):
    '''
    Calculate Gas transfer Velocity (cm/hr)
    Uses coefficients provided in Wanninkhof et al., (2014)
    
    Parameters
    ----------
    t: Int or Float 
        Temperature in Degrees Celcius

    Returns
    -------
    Sc : Float
        Gas Transfer Velocity in cm/hr
        Schmidt Coefficient
    '''
    #668 at 20C
    #Coefficients for Seawater CO2 Gas transfer rate
    a= 2116.8
    b=-136.25
    c= 4.7353
    d=-0.092307
    e=0.0007555
    Sc=a+(b*t)+(c*(t**2))+(d*(t**3)+(e*(t**4)))
    return Sc

def solubility(tk,s):
    '''
    Calculate CO2 solubility depending on temperature (Kelvin) and salinity
    Uses coefficients provided in Wanninkhof et al., (2014) and Weiss (1974)

    Parameters
    ----------
    tk : Int or Float
        Temperature in Kelvin (Degrees C +273.15)
    s : Int or Float
        Salinity in Practical Salinity Units (g/kg)        

    Returns
    -------
    ko : float
        Ko is the solubility of CO2 in seawater, given temperature and salinity.

    '''
    A1=-58.0931
    A2=90.5069
    A3=22.2940
    B1=0.027766
    B2=-0.025888
    B3=0.0050578

    ko=np.exp(A1+ A2 *(100/tk) + A3 *np.log(tk/100) + s*(B1 + B2 *(tk/100) + B3 *(tk/100)**2))

    return ko

def carbon_flux(s,t,u,pco2sw,pco2atm,dpco2=None,wsheight=None):
    '''
    Parameters
    ----------
    s : Int or Float
        Salinity in PSU (parts per thousand, g/kg).
    t :  Int or Float
        Temperature in Celcius
        (Will conver to kelvin automagically)
    u :  Int or Float
        Windspeed at 10m. (if wsheight!=None, will convert to 10m)
    pco2sw :  Int or Float
        pCO2 in SeaWater (ppm / uatm)
    pco2atm :  Int or Float
        pCO2 in Atmosphere (ppm / uatm)
    dpco2 :  Int or Float, optional
        Delta pCO2 (Difference, atm-sw; ppm/uatm)
        The default is None.
        if Not none, use instead of pco2sw and pco2atm
    wsheight : Int or Float
        Height of windspeed measurement (m)
        The default is None.
        if Not none, convert windspeed to U @ 10m
        using equation 1 in Sutton et al., (2017)
        
    Returns
    -------
    F : Float
        Co2 air-sea flux in mol / m2 / year
    Fa : Float
        Approximation: Co2 air-sea flux in mol/m2/year 
    '''

    tk=t+273.15 #Convert temperature to Kelvin
    
    if type(wsheight)==type(None):
        u10=u
    else:
        u10=(u) / 1 + ((np.sqrt(0.0011))/0.4)*np.log(wsheight/10) 
         #0.0011 cd10 drag coef and 0.4 von Karmans constant eq 1 Sutton 2017.
   
    if type(dpco2)==type(None):
        dpco2=pco2sw - pco2atm #Delta PCO2 between ocean and atmosphere. 
    
    Sc=schmidt(t) #Calculate Schmidt number from the functions above
    
    ko=solubility(tk,s)                     #Dimensionless solubility (Mol/k/atm)
    k=0.251*(u10**2)*(( Sc / 660)**-0.5)    #cm h –1
    #k1=0.251*(u10**2)*np.sqrt(660/Sc)      #Alternative calculation method as per oceanflux engine github.

    Fa=(7.7*(10**-4)*(u10**2)*dpco2)        #The approximation formula Wannikof 2014. 
    
    #Need to convert units into correct format as per comments here. 
    k_conv=(k*24*365)/100                    #cm/hr to m/yr              #Conversion = *24*365/100
    ko_conv=ko*1000                          #mol/atm/l to mol/atm/m3    #Conversion = *1000
    dpco2_conv=dpco2 /1000000                #uatm to atm                #Conversion = *1 000 000
    
    F=k_conv*ko_conv*dpco2_conv
    
    return F,Fa                           #mol/m2/year C
    

#Run tests
def schmidt_test():
    t=20
    ans=schmidt(t)
    assert int(schmidt(20))==668 #As per Sutton, Wannikof and Weiss.
    
def solubility_test():
    #The calues in Weiss 1974 were actually *10**-2
    assert np.round(solubility(tk=0+273.15,s=0),3)==np.round(7.758*10**-2,3) #Incorrect Ko solubility coeff
    assert np.round(solubility(tk=40+273.15,s=40),3)==np.round(2.044*10**-2,3)#Incorrect Ko solubility coeff
    assert np.round(solubility(tk=20+273.15,s=20),3)==np.round(3.562*10**-2,3) #Incorrect Ko solubility coeff

#The Schmidt and Solubulity functions work as expected.


def npp_to_moles(v):
    #mg C/m2/year
    #mol/m2/year #CO2
    molarmassC=12
    moles=(v/1000)/molarmassC #mol C / m2 / day
    return moles

def moles_to_carbon(moles):
    #In: moles Carbon/Co2
    #Out: g/C

    molarmassC=12
    grams=moles*molarmassC
    return grams


#s,t (c),u,pco2s,pco2a,dpco2
#print(carbon_flux(33,20,6,425,400))
#print(carbon_flux(35,10,10,475,400))

if __name__ == "__main__": 
    pass 



