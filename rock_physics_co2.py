#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 16:32:19 2023
Modified for Kimbelina 2 files
05/14/2025

Instructions:
python rock_physics_co2.py testdir mod1_090y.csv baseline_y0.csv

© 2025. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

 

"""
"""
Created on Tue Mar  9 17:11:14 2021
where material is formation (layer) ID
   Pre-Etchegoin: PE - sandstone
   Etchegoin: ET - sandstone?
   Macoma-Chanac: MC - shale
   Santa Margarita: SM  - sandstone
   Stevens: SV - sand
   Fruitvale-Round Mountain: RM - shale
   Olcese: OL -sand
   Temblor-Freeman: TF - CAP, shale
   Vedder: VC, VS - reservoir - sand
   Tumey-Eocence: TM - shale
   Last char: "C" = Clay and "F" = fault

@author: NEALA CREASY

Citations: 
    Wang, Z., Harbert, W. P., Dilmore, R. M., & Huang, L. (2018). 
    Modeling of time-lapse seismic monitoring using CO2 leakage simulations 
    for a model CO2 storage site with realistic geology: Application 
    in assessment of early leak-detection capabilities. 
    International Journal of Greenhouse Gas Control, 76, 39-52.
    
    Creasy, N., Huang, L., Gasperikova, E., Harbert, W., Bratton, T., & Zhou, Q. (2024). 
    CO2 rock physics modeling for reliable monitoring of geologic carbon storage. 
    Communications Earth & Environment, 5(1), 333.
"""

import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI
import math
import argparse

def run_vpvs(x,y,z,Pressure,P_initial,Tem,Saturation,Salt,porosity,DATADIR,filename):
    #print('good')
    #scalar values 
    rock = ['sandstone']
    #select the density of the dry rock matrix, sandstones and shales can vary
    #Input density of density of quartz (sandstone) 1.61-2.76, SHALE - 2.06-2.67    
    rho_dry = 2.6 #g/cm3 
    #Other input variables:
    v_qua = 0.7 #quartz/clay ratio, shale =< 0.5, sandstone ~ 0.7, so 0.7 = 70% quartz, 30% clay 
    pratio = 0.2 #Poisson ratio, can vary but shales and sandstones can be 0.2
    rho_matrix=rho_dry
    rho_c = rho_dry
    phi_cri=0.4  # Critical porosity of sandstone (from Nur et al. (1998)) 
    
    n = len(x)
    
    prefix = 'vel_pres_temp_sg_'
    fw2 = open(DATADIR +'/'+ prefix + filename, 'w')
    #col_values =['x','y','z','pres','temp','sg','salt','Vp','Vs','rhosat']
    fw2.write('index'+','+'x'+','+'y'+','+'z'+','+'pres'+','+'temp'+
                  ','+'sg'+','+'S'+','+'Vp'+','+'Vs'+','+'rho_sat' +','+
                  'porosity'+','+'Ksat'+','+'musat'+','+'znew'+','+'p_eff'+ '\n')
    
    for i in range(n):
        if (i % int(0.05*n)  == 0):
                print(int(math.ceil(i/n*100)),'%')
        
        T = Tem[i] # temperature C
        P = Pressure[i]
        P_i = P_initial[i]
        Sg = Saturation[i]
        S = Salt[i]
        poro = porosity[i]
        x_temp = x[i]
        y_temp = y[i]
        h = z[i]
        
        #1 Brine Batzle and Wang (1992) 
        # eq 27a:
        rho_water = 1+10**(-6)*(-80*Tem-3.3*Tem**2+0.00175*Tem**3 +
                        489*P-2*Tem*P+0.016*Tem**2*P-1.3*10**(-5)*Tem**3*P -
                        0.333*P**2-0.002*Tem*P**2)  # Density of pure water (g/cm3)

        #OUTPUT 1: rho_brine
        # eq 27b:
        rho_brine = rho_water+S*(0.668+0.44*S+1e-6*(300*P-2400*P*S +
                         Tem*(80+3*Tem-3300*S-13*P+47*P*S)))  # Density of brine (g/cm3)        

        # Table 1:
        wa = np.matrix([[1.40285e+03,  1.52400e+00,  3.43700e-03, -1.19700e-05],
                       [4.87100e+00, -1.11000e-02,  1.73900e-04, -1.62800e-06],
                        [-4.78300e-02,  2.74700e-04, -2.13500e-06,  1.23700e-08],
                        [1.48700e-04, -6.50300e-07, -1.45500e-08,  1.32700e-10],
                        [-2.19700e-07,  7.98700e-10,  5.23000e-11, -4.61400e-13]])
        wa_velocity = 0

        # Eq 28:
        for p in range(5):
            for j in range(4):
                wa_velocity += wa[p,j]*T**p*P**j  #P-wave velocity in pure water (m/s)
                
        #OUTPUT 2: brine_velocity
        # eq 29:
        brine_velocity = wa_velocity+S*(1170-9.6*Tem+0.055*Tem**2 -
                                8.5*10**(-5)*Tem**3+2.6*P-0.0029*Tem*P -
                                0.0476*P**2)+S**1.5*(780-10*P+0.16*P**2)-1820*S**2  # P-wave velocity of brine (m/s)
        k_brine=rho_brine*brine_velocity**2*10**(-6) ###Bulk modulus of brine (GPa)      
        
        
        # 2 CONFINING AND DIFFERENTIAL PRESSURE #########################
        # to find confining pressure, we assume hydrostatic pressure based on 
        # the minimum pore pressure and assume fresh water
        Pmin = np.min(P_initial)*1e6 #pascal
        dens_ratio = 1000/(rho_c*1000) # density ratio for top of model
        Pc_0 = Pmin/dens_ratio # minimum confining pressure based on hydrostatic load
        g=9.81 #m/s2
        h_0 = Pmin/(g*1000)
        h_new = h+h_0 #new adjust depth
        #NEED DENSITY, PORE PRESSURE
        #get average density of fluid in model based on pressure gradient of fluid
        #pressure gradient, pascal/meter
        grad4 = (np.max(P_initial)*1e6-np.min(P_initial)*1e6)/(np.max(z)-np.min(z))
        ave_rh = grad4/g #average density of entire liquid in model
        pc_model = P_i/(ave_rh/(rho_c*1000))*10**6 #pascal, uses initial pore pressure to find confining pressure
        pc = Pc_0+pc_model #confinig pressure, pascals
        P_Pa=P*10**6 #convert pore pressure to pascal from MPa
        pd_mpa = (pc-P_Pa)*10**-6 #differential in MPA
        pd1=(pc-P_Pa)*10**-9 # differential pressure (GPa), pd = pc - p
        #print(P,h_0,h_new,h)
        
        #CALCULATE MINERAL PROPERTIES: Voigt-Reuss-Hill averaging (Hill, 1952) 
        k_quartz = 37.0 #Bulk modulus of quartz mineral (GPa), Kimizuka et al., 2007
        k_clay = 21.0
        v_clay = 1 - v_qua
        mu_quartz=44.0 #Shear modulus of quartz mineral (GPa)line 38
        mu_clay=7.0 #Shear modulus of feldspar mineral (GPa)
        mu_up=np.amax([mu_quartz,mu_clay])
        mu_lo=np.amin([mu_quartz,mu_clay])
        k_up=np.amax([k_quartz,k_clay])
        k_lo=np.amin([k_quartz,k_clay])
        k_HSlo=(v_qua/(k_quartz+4.0/3*mu_lo)+v_clay/(k_clay+4.0/3*mu_lo))**-1-4.0/3*mu_lo  ##GPa
        k_HSup=(v_qua/(k_quartz+4.0/3*mu_up)+v_clay/(k_clay+4.0/3*mu_up))**-1-4.0/3*mu_up  ##GPa
        ksi_up=mu_up/6*((9*k_up+8*mu_up)/(k_up+2*mu_up))
        ksi_lo=mu_lo/6*((9*k_lo+8*mu_lo)/(k_lo+2*mu_lo))
        mu_HSup=(v_qua/(mu_quartz+ksi_up)+v_clay/(mu_clay+ksi_up))**(-1)-ksi_up  ##GPa
        mu_HSlo=(v_qua/(mu_quartz+ksi_lo)+v_clay/(mu_clay+ksi_lo))**(-1)-ksi_lo  ##GPa
        k_grains=(k_HSup+k_HSlo)/2  ## Bulk modulus of the grains (Gpa)
        mu_grains=(mu_HSup+mu_HSlo)/2  ## Shear modulus of the grains (Gpa)
        
        
        sporo = poro
        phic0 = 5.073
        thetac = 1607
        thetacu = 1339
        tCd = 0.085
        cpporo=poro+phic0*math.exp(-tCd*pd_mpa)/100
        
        phic0co2 = 8.15
        tCdco2 = 0.08
        thetacc = 1193
        thetacuc = 824
        Dp,D = 0.085,0.08 #pre,post
        kdiff = 0.1147 #difference in KDRS
        udiff = 0.2544
        cctu = 3.174/1000 #(C-C)*thetasu
        cct = 1.0929/1000
        
        if Sg>0:
            cporo = poro + phic0co2*math.exp(-tCdco2*pd_mpa)/100
            poro = cporo
            
        else:
            poro = cpporo
        
            
        
        ###### 10) CO2 Supercritical Properties (Friedman, 1963; Carcione et al., 2006)
        #R=8.314 #gas constant J/(mol*K)
        Tem_K=T+273.15 #Kelvin
        vel_co2=PropsSI('speed_of_sound', 'T', Tem_K, 'P', P_Pa, 'CO2') # m/s
        rho_gas=PropsSI('D', 'T', Tem_K, 'P', P_Pa, 'CO2')*10**-3 # g/cm3
        k_co2=rho_gas*vel_co2**2*10**(-6) #CO2 bulk modulus in GPa
        Sw=1-Sg # Water(Brine) saturation
        K_fl=1/(Sw/k_brine+Sg/k_co2) # Bulk modulus of the fluid phase (GPa)
        rho_fl=Sw*rho_brine+Sg*rho_gas # Density of the fluid phase (g/cm3)
        
        
        
        ###### 11) Rho SAT   ######################################################
        rho_sat=poro*rho_fl+(1-poro)*rho_matrix # Density of rock after fluid substitution (g/cm3)
        
        ###### 12b) Hashin Strikman + CP ###################################################
        C=2.8/phi_cri  # from Murphy (1982)
        Khm=(C**2*(1-phi_cri)**2*mu_grains**2*pd1/(18*np.pi**2*(1-pratio)**2))**(1.0/3) 
        # Bulk modulus at the critical porosity (GPa)
        Ghm=(5-4.0*pratio)/(5*(2-pratio))*(3*C**2*(1-phi_cri)**2*
                                           mu_grains**2*pd1/
                                           (2*np.pi**2*(1-pratio)**2))**(1.0/3) 
        # Shear modulus at the critical porosity (GPa)
        Rphi=sporo/phi_cri  #must use stiff porosity
        K_frame=(Rphi/(Khm+4.0/3*Ghm)+(1-Rphi)/(k_grains+4.0/3*Ghm))**(-1)-4.0/3*Ghm
        mu_frame=(Rphi/(Ghm+Ghm/6.0*(9*Khm+8*Ghm)/(Khm+2*Ghm))+
                  (1-Rphi)/(mu_grains+Ghm/6*(9*Khm+8*Ghm)/(Khm+2*Ghm)))**(-1)-Ghm/6*(9*Khm+8*Ghm)/(Khm+2*Ghm)
            
        
        #K_frame = K_frame*(1-thetac*phic0*1e-4*math.exp(-tCd*pd_mpa))
        
        if Sg > 0:
            #k_change = (kpost-kpre)/kpre
            #K_frame = K_frame*(1+k_change)
            K_frame = K_frame*(1-kdiff)*(1+cct*pd_mpa-thetacc*phic0co2*1e-4*math.exp(-tCdco2*pd_mpa))
            mu_frame = mu_frame*(1-udiff)*(1+cctu*pd_mpa-thetacuc*phic0co2*1e-4*math.exp(-tCdco2*pd_mpa))
        else:
            K_frame = K_frame*(1-thetac*phic0*1e-4*math.exp(-tCd*pd_mpa))
            mu_frame = mu_frame*(1-thetacu*phic0*1e-4*math.exp(-tCd*pd_mpa))
                    
        
        # - Gassmann
        K_sat=K_frame+(1-K_frame/k_grains)**2/(poro/K_fl+(1-poro)/k_grains-K_frame/(k_grains**2))  ### Bulk modulus of rock after fluid substitution (Gpa)
        
        mu_sat = mu_frame
        
        ###### 14) Find Vp and Vs ##############################################
        #print "Bulk modulus of rock after fluid substitution =", K_sat, "Gpa"
        pVelocity=((K_sat*10**9+4.0/3*mu_sat*10**9)/(rho_sat*10**3))**0.5  ### P-wave velocity (m/s)
        sVelocity=((mu_sat*10**9)/(rho_sat*10**3))**0.5  ### S-wave velocity (m/s)
        
        rho_sat_out = rho_sat
        fw2.write(str(i)+',' 
                      + format(x[i],'.2f')+',' + format(y[i],'.2f')
                      +','+ format(z[i],'.2f')+',' + format(P,'.2f')
                      +',' + format(T,'.2f')+',' + format(Sg,'.4f')
                      +',' + format(S,'.9f')+','
                      + format(pVelocity, '.2f')
                      +',' + format(sVelocity, '.2f')+','+ format(rho_sat_out, '.4f')
                      +','+ format(poro,'.4f') 
                      +',' + format(K_sat,'.4f') +',' + format(mu_sat,'.4f') 
                      +',' + format(h_new,'.2f')+','+format(pd_mpa,'.2f')+'\n')
        
    fw2.close()
  

        
    

def read_file(DATADIR,filename,baseline):
    ######## LOADS DATA FILES
    data_list=DATADIR+'/'+filename
    base_list=DATADIR+'/'+baseline
    print(data_list)
    print(base_list)
    base = pd.read_csv(base_list, sep=',')
    #print(base)
    data1 = pd.read_csv(data_list, sep=',')
    #print(data1)
    #check length
    count = len(base)
    print('Number of Baseline data file lines to be processed', count)
    count2 = len(data1)
    print('Number of Simulation file lines to be processed', count2)
    #confirm that baseline and simulation files are same length
    if count != count2: 
           print('baseline and simulation files are different lengths')
           
    Tem = data1.temp
    P = data1.pres/(1.00000e+06) #PORE PRESSURE! MPA
    P_i = base.P/(1.00000e+06) #Baseline PORE PRESSURE! MPA
    S = (data1.salt) #salt content
    Sg = data1.sg #Gas concentration
    # Depth below surface is (z value) Depth units are in meters (double check)
    h = abs(data1.z - np.max(data1.z))  
    
    
    poro = base.porosity
    x = base.x
    y = base.y
    
    
    
    run_vpvs(x,y,h,P,P_i,Tem,Sg,S,poro,DATADIR,filename)
    
    


if __name__ == "__main__":
    # Process input arguments
    myname="rock_physics_co2"
    s="Mon Oct 23 2023"
    lcdate="(last change: "+s+")" 
    s="""Original version description: Caclculate seismic velocity model Vp and Vs from 
    the TOUGH reservoir model grid. 
    Save the model in a (zipped) csv file, one file per time step. 
    
    """
    tlog=["{} is a Python3 script for rock physics modeling of CSV-formatted data input file(s). ".format(myname)]
    tlog.append(s)
    s=" Authors: Neala Creasy (2023), LANL "
    tlog.append(s)
    
    parser = argparse.ArgumentParser(prog=myname,description=tlog[0]+tlog[1]+tlog[2]+lcdate,epilog=s)
    # TOUGH2 data file in following CSV format (example):
    #elem,x,y,z,pres,temp,sg,xco2\n 
    #AI01a 1,     0.00,     0.00,   329.75, 4.049016e+06, 30.83, 0.000, 0.0E+00\n 
    #BI01a 1,     0.00,     0.00,  1125.32, 1.182724e+07, 52.15, 0.000, 0.0E+00\n 
    # Main arguments
    s="Root pathname of data input file directories"
    parser.add_argument('DATADIR',help=s)
    s="File name"
    parser.add_argument('filename',help=s)
    s="Basline file name"
    parser.add_argument('baseline',help=s)
    #s="List of comma-separated CSV-formatted data input data files, data files have one header line"
    #parser.add_argument('DATAFILE',help=s)
 
    args = parser.parse_args() # parse all arguments
    #
    print("{} - Calculation of geophysical properties".format(myname))
    s="""Updated version {} combines seismic property and electrical-resistivity
    property generation {}.""".format(myname,lcdate)
    tlog.append(s)
    
    #read in files
    read_file(args.DATADIR,args.filename,args.baseline)
        

 
