# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 12:30:47 2016


Functions for the analysis of WRF microphysical pathways
@author: heikenfeld
"""

from collections import defaultdict


#Thompson Microphysics:
#    
#    !+---+-----------------------------------------------------------------+
#!.. Source/sink terms.  First 2 chars: "pr" represents source/sink of
#!.. mass while "pn" represents source/sink of number.  Next char is one
#!.. of "v" for water vapor, "r" for rain, "i" for cloud ice, "w" for
#!.. cloud water, "s" for snow, and "g" for graupel.  Next chars
#!.. represent processes: "de" for sublimation/deposition, "ev" for
#!.. evaporation, "fz" for freezing, "ml" for melting, "au" for
#!.. autoconversion, "nu" for ice nucleation, "hm" for Hallet/Mossop
#!.. secondary ice production, and "c" for collection followed by the
#!.. character for the species being collected.  ALL of these terms are
#!.. positive (except for deposition/sublimation terms which can switch
#!.. signs based on super/subsaturation) and are treated as negatives
#!.. where necessary in the tendency equations.
#!+---+-----------------------------------------------------------------+

#       
#    
List_Processes_Thompson_Mass=[
         'PRW_VCD',   #  Vapor->Water  
         'PRV_REV',   #  Vapor->Water  
         'PRR_WAU',   #  Vapor->Water  
         'PRR_RCW',   #  Vapor->Water  
         'PRR_RCS',   #  Vapor->Water  
         'PRR_RCG',   #  Rain->Graupel  
         'PRR_GML',   #  Vapor->Water  
         'PRR_RCI',   #  Vapor->Water          
         'PRI_INU',   #  Vapor->Water     
         'PRI_IHM',   #  Vapor->Water     
         'PRI_WFZ',   #  Vapor->Water    
         'PRI_RFZ',   #  Vapor->Water  
         'PRI_IDE',   #  Vapor->Water  
         'PRI_RCI',   #  Vapor->Water  
         'PRI_IHA',   #  Vapor->Water  
         'PRS_IAU',   #  Vapor->Water  
         'PRS_SCI',   #  Vapor->Water  
         'PRS_RCS',   #  Vapor->Water  
         'PRS_SCW',   #  Vapor->Water  
         'PRS_SDE',   #  Vapor->Water  
         'PRS_IHM',   #  Vapor->Water  
         'PRS_IDE',   #  Vapor->Water  
         'PRG_SCW',   #  Vapor->Water  
         'PRG_RFZ',   #  Vapor->Water  
         'PRG_GDE',   #  Vapor->Water  
         'PRG_GCW',   #  Vapor->Water  
         'PRG_RCI',   #  Vapor->Water  
         'PRG_RCS',   #  Vapor->Water  
         'PRG_RCG',   #  Rain->Graupel  
         'PRG_IHM'    #  Graupel->Ice                               
          ]    
          
Processes_Thompson_Mass_colors={}
Processes_Thompson_Mass_colors['PRW_VCD']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRV_REV']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRR_WAU']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRR_RCW']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRR_RCS']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRR_RCG']='slategrey'   #  Rain->Graupel  
Processes_Thompson_Mass_colors['PRR_GML']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRR_RCI']='slategrey'   #  Vapor->Water
Processes_Thompson_Mass_colors['PRI_INU']='slategrey'   #  Vapor->Water     
Processes_Thompson_Mass_colors['PRI_IHM']='slategrey'   #  Vapor->Water     
Processes_Thompson_Mass_colors['PRI_WFZ']='slategrey'   #  Vapor->Water    
Processes_Thompson_Mass_colors['PRI_RFZ']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRI_IDE']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRI_RCI']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRI_IHA']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRS_IAU']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRS_SCI']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRS_RCS']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRS_SCW']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRS_SDE']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRS_IHM']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRS_IDE']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRG_SCW']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRG_RFZ']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRG_GDE']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRG_GCW']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRG_RCI']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRG_RCS']='slategrey'   #  Vapor->Water  
Processes_Thompson_Mass_colors['PRG_RCG']='slategrey'   #  Rain->Graupel  
Processes_Thompson_Mass_colors['PRG_IHM']='slategrey'   #  Graupel->Ice                               
                    
  


List_Processes_Thompson_Number=[
         'PNC_WCD',
         'PNC_WAU',
         'PNC_RCW',
         'PNC_SCW',
         'PNC_GCW',
         'PNR_WAU',
         'PNR_RCS',
         'PNR_RCG',
         'PNR_RCI',
         'PNR_SML',
         'PNR_GML',
         'PNR_REV',
         'PNR_RCR',
         'PNR_RFZ',
         'PNI_INU',
         'PNI_IHM',
         'PNI_WFZ',
         'PNI_RFZ',
         'PNI_IDE',
         'PNI_RCI',
         'PNI_SCI',
         'PNI_IAU',
         'PNI_IHA',
         'PNA_RCA',
         'PNA_SCA',
         'PNA_GCA',      
         'PNA_RCA',
         'PNA_SCA',
         'PNA_GCA',             
         'PND_RCD',
         'PND_SCD',
         'PND_GCD'               
                              ]    

def load_wrf_thom_mass_proc(filename,add_coordinates=None,quantity='volume'):
    from wrfload import loadwrfcube, derivewrfcube
    Dict={}
    if quantity=='volume':
        rho=derivewrfcube(filename,'density')
        
    for i_process,process in enumerate(List_Processes_Thompson_Mass):
        #print(process)
        if (i_process==0):
            cube=loadwrfcube(filename,process,add_coordinates=None)
            cube.rename(process)
            #Cubelist.append(cube)
            Dict[process]=cube
            if add_coordinates=='pz':
               z_coord=cube.coord('geopotential')
               p_coord=cube.coord('pressure')
        else:
            cube=loadwrfcube(filename,process)
            cube.rename(process)
            if add_coordinates=='pz':
               cube.add_aux_coord(z_coord,(0,1,2,3))
               cube.add_aux_coord(p_coord,(0,1,2,3))
        #Cubelist.append(cube)
        if quantity=='volume':
            cube=cube*rho
        Dict[process]=cube
    #return Cubelist
    return Dict
 
def load_wrf_thom_number_proc(filename,add_coordinates=None,quantity='volume'):
    from wrfload import loadwrfcube, derivewrfcube
    Dict={}
    if quantity=='volume':
        rho=derivewrfcube(filename,'density')

    for i_process,process in enumerate(List_Processes_Thompson_Number):
        #print(process)
        if (i_process==0):
            cube=loadwrfcube(filename,process,add_coordinates=None)
            cube.rename(process)
            #Cubelist.append(cube)
            Dict[process]=cube
            if add_coordinates=='pz':
               z_coord=cube.coord('geopotential')
               p_coord=cube.coord('pressure')
        else:
            cube=loadwrfcube(filename,process)
            cube.rename(process)
            if add_coordinates=='pz':
                cube.add_aux_coord(z_coord,(0,1,2,3))
                cube.add_aux_coord(p_coord,(0,1,2,3))
        #Cubelist.append(cube)
        if quantity=='volume':
            cube=cube*rho

        Dict[process]=cube
    #return Cubelist
    return Dict                     
                              
#
#! MICROPHYSICAL PROCESSES
#
#     REAL, DIMENSION(KTS:KTE) ::  NSUBC     ! LOSS OF NC DURING EVAP
#     REAL, DIMENSION(KTS:KTE) ::  NSUBI     ! LOSS OF NI DURING SUB.
#     REAL, DIMENSION(KTS:KTE) ::  NSUBS     ! LOSS OF NS DURING SUB.
#     REAL, DIMENSION(KTS:KTE) ::  NSUBR     ! LOSS OF NR DURING EVAP
     #     REAL, DIMENSION(KTS:KTE) ::  PRD       ! DEP CLOUD ICE
     #     REAL, DIMENSION(KTS:KTE) ::  PRE       ! EVAP OF RAIN
     #     REAL, DIMENSION(KTS:KTE) ::  PRDS      ! DEP SNOW
#     REAL, DIMENSION(KTS:KTE) ::  NNUCCC    ! CHANGE N DUE TO CONTACT FREEZ DROPLETS
     #     REAL, DIMENSION(KTS:KTE) ::  MNUCCC    ! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
     #     REAL, DIMENSION(KTS:KTE) ::  PRA       ! ACCRETION DROPLETS BY RAIN
     #     REAL, DIMENSION(KTS:KTE) ::  PRC       ! AUTOCONVERSION DROPLETS
    #     REAL, DIMENSION(KTS:KTE) ::  PCC       ! COND/EVAP DROPLETS
#     REAL, DIMENSION(KTS:KTE) ::  NNUCCD    ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
    #     REAL, DIMENSION(KTS:KTE) ::  MNUCCD    ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
    #     REAL, DIMENSION(KTS:KTE) ::  MNUCCR    ! CHANGE Q DUE TO CONTACT FREEZ RAIN
#     REAL, DIMENSION(KTS:KTE) ::  NNUCCR    ! CHANGE N DUE TO CONTACT FREEZ RAIN
#     REAL, DIMENSION(KTS:KTE) ::  NPRA      ! CHANGE IN N DUE TO DROPLET ACC BY RAIN
#     REAL, DIMENSION(KTS:KTE) ::  NRAGG     ! SELF-COLLECTION/BREAKUP OF RAIN
#     REAL, DIMENSION(KTS:KTE) ::  NSAGG     ! SELF-COLLECTION OF SNOW
#     REAL, DIMENSION(KTS:KTE) ::  NPRC      ! CHANGE NC AUTOCONVERSION DROPLETS
#     REAL, DIMENSION(KTS:KTE) ::  NPRC1      ! CHANGE NR AUTOCONVERSION DROPLETS
    #     REAL, DIMENSION(KTS:KTE) ::  PRAI      ! CHANGE Q ACCRETION CLOUD ICE BY SNOW
    #     REAL, DIMENSION(KTS:KTE) ::  PRCI      ! CHANGE Q AUTOCONVERSIN CLOUD ICE TO SNOW
    #     REAL, DIMENSION(KTS:KTE) ::  PSACWS    ! CHANGE Q DROPLET ACCRETION BY SNOW
#     REAL, DIMENSION(KTS:KTE) ::  NPSACWS   ! CHANGE N DROPLET ACCRETION BY SNOW
    #     REAL, DIMENSION(KTS:KTE) ::  PSACWI    ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
#     REAL, DIMENSION(KTS:KTE) ::  NPSACWI   ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
#     REAL, DIMENSION(KTS:KTE) ::  NPRCI     ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
#     REAL, DIMENSION(KTS:KTE) ::  NPRAI     ! CHANGE N ACCRETION CLOUD ICE
#     REAL, DIMENSION(KTS:KTE) ::  NMULTS    ! ICE MULT DUE TO RIMING DROPLETS BY SNOW
#     REAL, DIMENSION(KTS:KTE) ::  NMULTR    ! ICE MULT DUE TO RIMING RAIN BY SNOW
    #     REAL, DIMENSION(KTS:KTE) ::  QMULTS    ! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
    #     REAL, DIMENSION(KTS:KTE) ::  QMULTR    ! CHANGE Q DUE TO ICE RAIN/SNOW
#     REAL, DIMENSION(KTS:KTE) ::  PRACS     ! CHANGE Q RAIN-SNOW COLLECTION                          to ICE or GRAUP ??????
#     REAL, DIMENSION(KTS:KTE) ::  NPRACS    ! CHANGE N RAIN-SNOW COLLECTION
    #     REAL, DIMENSION(KTS:KTE) ::  PCCN      ! CHANGE Q DROPLET ACTIVATION
    #     REAL, DIMENSION(KTS:KTE) ::  PSMLT     ! CHANGE Q MELTING SNOW TO RAIN
    #     REAL, DIMENSION(KTS:KTE) ::  EVPMS     ! CHNAGE Q MELTING SNOW EVAPORATING
#     REAL, DIMENSION(KTS:KTE) ::  NSMLTS    ! CHANGE N MELTING SNOW
#     REAL, DIMENSION(KTS:KTE) ::  NSMLTR    ! CHANGE N MELTING SNOW TO RAIN
#! HM ADDED 12/13/06
    #     REAL, DIMENSION(KTS:KTE) ::  PIACR     ! CHANGE QR, ICE-RAIN COLLECTION
#     REAL, DIMENSION(KTS:KTE) ::  NIACR     ! CHANGE N, ICE-RAIN COLLECTION
    #     REAL, DIMENSION(KTS:KTE) ::  PRACI     ! CHANGE QI, ICE-RAIN COLLECTION
    #     REAL, DIMENSION(KTS:KTE) ::  PIACRS     ! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
#     REAL, DIMENSION(KTS:KTE) ::  NIACRS     ! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
    #     REAL, DIMENSION(KTS:KTE) ::  PRACIS     ! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
    #     REAL, DIMENSION(KTS:KTE) ::  EPRD      ! SUBLIMATION CLOUD ICE
    #     REAL, DIMENSION(KTS:KTE) ::  EPRDS     ! SUBLIMATION SNOW
#! HM ADDED GRAUPEL PROCESSES
    #     REAL, DIMENSION(KTS:KTE) ::  PRACG    ! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
    #     REAL, DIMENSION(KTS:KTE) ::  PSACWG    ! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
    #     REAL, DIMENSION(KTS:KTE) ::  PGSACW    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
    #     REAL, DIMENSION(KTS:KTE) ::  PGRACS    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
    #     REAL, DIMENSION(KTS:KTE) ::  PRDG    ! DEP OF GRAUPEL
    #     REAL, DIMENSION(KTS:KTE) ::  EPRDG    ! SUB OF GRAUPEL
    #     REAL, DIMENSION(KTS:KTE) ::  EVPMG    ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
    #     REAL, DIMENSION(KTS:KTE) ::  PGMLT    ! CHANGE Q MELTING OF GRAUPEL
#     REAL, DIMENSION(KTS:KTE) ::  NPRACG    ! CHANGE N COLLECTION RAIN BY GRAUPEL
#     REAL, DIMENSION(KTS:KTE) ::  NPSACWG    ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
#     REAL, DIMENSION(KTS:KTE) ::  NSCNG    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
#     REAL, DIMENSION(KTS:KTE) ::  NGRACS    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
#     REAL, DIMENSION(KTS:KTE) ::  NGMLTG    ! CHANGE N MELTING GRAUPEL
#     REAL, DIMENSION(KTS:KTE) ::  NGMLTR    ! CHANGE N MELTING GRAUPEL TO RAIN
#     REAL, DIMENSION(KTS:KTE) ::  NSUBG    ! CHANGE N SUB/DEP OF GRAUPEL
#     REAL, DIMENSION(KTS:KTE) ::  PSACR    ! CONVERSION DUE TO COLL OF SNOW BY RAIN
#     REAL, DIMENSION(KTS:KTE) ::  NMULTG    ! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
#     REAL, DIMENSION(KTS:KTE) ::  NMULTRG    ! ICE MULT DUE TO ACC RAIN BY GRAUPEL
    #     REAL, DIMENSION(KTS:KTE) ::  QMULTG    ! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
    #     REAL, DIMENSION(KTS:KTE) ::  QMULTRG    ! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL




#Morrison microphysics

#Morrison microphysics
Proclist_Morr_mass_load=[
    'PRD',
    'PRE',
    'PRDS',
    'PRA',
    'PRC',
    'PCC',
    #'PCCN', currently wrongly written out by the pathway analysis
    'PSMLT',
    'EVPMS',
    'QMULTS',
    'QMULTR',
    'PRACS',
    'PSACWG',
    'PGSACW',
    'PGRACS',
    'PRDG',
    'EPRDG',
    'EVPMG',
    'PGMLT',
    'PRACI',
    'PIACRS',
    'PRACIS',
    'EPRD',
    'EPRDS',
    'PRACG',
    'QMULTG',
    'MNUCCR',
    'QMULTRG',
    'PRAI',
    'PRCI',
    'PSACWS',
    'PIACR',
    'PSACWI',
    'PSACR']

Proclist_Morr_mass=list(Proclist_Morr_mass_load)
Proclist_Morr_mass.append('EPCC')

Morr_Processes_colors=defaultdict(dict)
Morr_Processes_colors['PRD']='lightseagreen'
Morr_Processes_colors['PRE']='red'
Morr_Processes_colors['PRDS']='darkorange'
Morr_Processes_colors['PRA']='darkred'
Morr_Processes_colors['PRC']='orange'
Morr_Processes_colors['PCC']='magenta'
#Morr_Processes_colors['PCCN']='blue'
Morr_Processes_colors['PSMLT']='darkblue'
Morr_Processes_colors['EVPMS']='azure'
Morr_Processes_colors['QMULTS']='cyan'
Morr_Processes_colors['QMULTR']='green'
Morr_Processes_colors['PRACS']='darkgreen'
Morr_Processes_colors['PSACWG']='palegreen'
Morr_Processes_colors['PGSACW']='gray'
Morr_Processes_colors['PGRACS']='black'
Morr_Processes_colors['PRDG']='springgreen'
Morr_Processes_colors['EPRDG']='coral'
Morr_Processes_colors['EVPMG']='sage'
Morr_Processes_colors['PGMLT']='mediumpurple'
Morr_Processes_colors['PRACI']='lightsteelblue'
Morr_Processes_colors['PIACRS']='darkslategrey'
Morr_Processes_colors['PRACIS']='skyblue'
Morr_Processes_colors['EPRD']='beige'
Morr_Processes_colors['EPRDS']='saddlebrown'
Morr_Processes_colors['PRACG']='violet'
Morr_Processes_colors['QMULTG']='pink'
Morr_Processes_colors['QMULTRG']='indigo'
Morr_Processes_colors['MNUCCR']='lightcyan'
Morr_Processes_colors['PRAI']='lime'
Morr_Processes_colors['PRCI']='peru'
Morr_Processes_colors['PSACWS']='maroon'
Morr_Processes_colors['PIACR']='black'
Morr_Processes_colors['PSACWI']='gold'
Morr_Processes_colors['PSACR']='lightgray'
Morr_Processes_colors['EPCC']='lightblue'

def load_wrf_morr_mass_proc(filename,add_coordinates=None,quantity='mixing ratio',slice_time=slice(None)):
    import numpy as np
    from wrfload import loadwrfcube, derivewrfcube
    Proclist=Proclist_Morr_mass_load
    Dict={}
    print('start loading processes')

    if quantity=='volume':
        print('start calculating density')
        rho=derivewrfcube(filename,'density',slice_time=slice_time)
    for i_process,process in enumerate(Proclist):
        print(process)
        if (i_process==0):
            cube=loadwrfcube(filename,process+'3D',add_coordinates=add_coordinates,slice_time=slice_time)
            #print(cube)
            #print(add_coordinates)
            if add_coordinates=='pz':
                z_coord=cube.coord('geopotential')
                p_coord=cube.coord('pressure')
            if quantity=='volume':
                cube=cube*rho
            cube.rename(process)

            Dict[process]=cube

        else:
            #print(process)
            if process=='PCC':
                cube=loadwrfcube(filename,process+'3D',slice_time=slice_time)
                cube1=cube[:]
                cube2=cube[:]
                cube1.data=np.clip(cube1.data,a_min=-np.inf,a_max=0)            
                if add_coordinates=='pz':
                    cube1.add_aux_coord(z_coord,(0,1,2,3))
                    cube1.add_aux_coord(p_coord,(0,1,2,3))
                if quantity=='volume':
                    cube1=cube1*rho                
                    
                cube1.rename('PCC')
                Dict['PCC']=cube1
                cube2=cube[:]
                cube2.data=np.clip(cube1.data,a_min=0,a_max=np.inf)
                if add_coordinates=='pz':
                    cube2.add_aux_coord(z_coord,(0,1,2,3))
                    cube2.add_aux_coord(p_coord,(0,1,2,3))
                if quantity=='volume':
                    cube2=cube2*rho
                cube2.rename('EPCC')
                Dict['EPCC']=cube2
    
            else:
                cube=loadwrfcube(filename,process+'3D',slice_time=slice_time)
                if add_coordinates=='pz':
                    cube.add_aux_coord(z_coord,(0,1,2,3))
                    cube.add_aux_coord(p_coord,(0,1,2,3))
                #Cubelist.append(cube)
                if quantity=='volume':
                    cube=cube*rho
                cube.rename(process)
                print(cube)
                Dict[process]=cube

    #return Cubelist
    return Dict
    
def load_wrf_morr_mass_proc_individual(filename,process,add_coordinates=None,quantity='mixing ratio',slice_time=slice(None)):
    import numpy as np
    from wrfload import loadwrfcube, derivewrfcube
    Dict={}
    if quantity=='volume':
        rho=derivewrfcube(filename,'density',slice_time=slice_time)

    if process=='PCC':
        cube=loadwrfcube(filename,process+'3D',add_coordinates=add_coordinates,slice_time=slice_time)
        cube1=cube.copy()
        cube2=cube.copy()
        cube1.data=np.clip(cube1.data,a_min=-np.inf,a_max=0)            
        if quantity=='volume':
            cube1=cube1*rho                
        cube1.rename('PCC')
        Dict['PCC']=cube1
        cube2=cube[:]
        cube2.data=np.clip(cube.data,a_min=0,a_max=np.inf)
        if quantity=='volume':
            cube2=cube2*rho
        cube2.rename('EPCC')
        Dict['EPCC']=cube2

    else:
        cube=loadwrfcube(filename,process+'3D',add_coordinates=add_coordinates,slice_time=slice_time)
        if add_coordinates=='pz':
            cube.add_aux_coord(z_coord,(0,1,2,3))
            cube.add_aux_coord(p_coord,(0,1,2,3))
        #Cubelist.append(cube)
        if quantity=='volume':
            cube=cube*rho
        cube.rename(process)

        Dict[process]=cube

    #return Cubelist
    return Dict    
    
    

Proclist_Morr_number=['NSUBC',
    'NSUBI',
    'NSUBS',
    'NSUBR',
    'NNUCCC',
    'MNUCCC',
    'NNUCCD',
    'MNUCCD',
    'NNUCCR',
    'NPRA',
    'NRAGG',
    'NSAGG',
    'NPRC',
    'NPRC1',
    'NPSACWS',
    'NPSACWI',
    'NPRCI',
    'NPRAI',
    'NMULTS',
    'NMULTR',
    'NPRACS',
    'NSMLTR',
    'NSMLTS',
    'NIACR',
    'NPRACG',
    #'NPSACWG',
    'NSCNG',
    'NGRACS',
    'NGMLTG',
    'NGMLTR',
    'NSUBG',
    'NMULTG',
    'NMULTRG,'    
    'NIACRS']
    
def load_wrf_morr_num_proc(filename,add_coordinates=None,quantity='volume'):
    from wrfload import loadwrfcube, derivewrfcube
    Dict={}
    if quantity=='volume':
        rho=derivewrfcube(filename,'density')

 
    #Cubelist=[]
    Dict={}
    for i_process,process in enumerate(Proclist_Morr_number):
        if (i_process==0):
            cube=loadwrfcube(filename,process+'3D',add_coordinates=add_coordinates)
            cube.rename(process)
            #Cubelist.append(cube)
            if add_coordinates=='pz':
                z_coord=cube.coord('geopotential')
                p_coord=cube.coord('pressure')
        else:
            if process=='NSMLTR':
                cube=loadwrfcube(filename,process)
            else:
                cube=loadwrfcube(filename,process+'3D')
            cube.rename(process)
            if add_coordinates=='pz':
                cube.add_aux_coord(z_coord,(0,1,2,3))
                cube.add_aux_coord(p_coord,(0,1,2,3))

        #Cubelist.append(cube)
        if quantity=='volume':
            cube=cube*rho

        Dict[process]=cube
    #return Cubelist
    return Dict

Hydropath_list=[
    'VAPORCLOUD',
    'VAPORRAIN',
    'VAPORICE',
    'VAPORSNOW',
    'VAPORGRAUP',
    'CLOUDRAIN',
    'CLOUDICE',
    'CLOUDSNOW',
    'CLOUDGRAUP',
    'RAINICE',
    'RAINSNOW',
    'RAINGRAUP',
    'ICESNOW',
    'ICEGRAUP',
    'SNOWGRAUP']    
    
def calculate_wrf_morr_path_hydrometeors(filename,slice_time=slice(None)):
    Dict={}
    #Cubelist=[]
    for path in Hydropath_list:
        cube=calculate_wrf_morr_path(filename,path,slice_time=slice_time)
        #Cubelist.append(cube)
        Dict[path]=cube
    #return Cubelist
    return Dict
    
Phasepath_list=['vaporliquid','vaporfrozen','liquidfrozen']

def calculate_wrf_morr_path_phases(filename,slice_time=slice(None)):
    Dict={}
    #Cubelist=[]
    for path in Phasepath_list:
        print('loading ',  path)
        cube=calculate_wrf_morr_path(filename,path,slice_time=slice_time)
        print(path, ' loaded')

        #Cubelist.append(cube)
        Dict[path]=cube
        print(Dict[path].data)
    #return Cubelist
    return Dict
    
def calculate_wrf_thompson_path(filename,path,add_coordinates=None,quantity='volume'):
    if (path=='processes_mass'):
        out=load_wrf_thom_mass_proc(filename,add_coordinates)
    if (path=='processes_number'):
        out=load_wrf_thom_number_proc(filename,add_coordinates)
    else:
        print('option not avaliable')
    return out   

        
def calculate_wrf_morr_path(filename,path,add_coordinates=None,quantity='volume',slice_time=slice(None)):
    if (path=='processes_mass'):
        out=load_wrf_morr_mass_proc(filename,add_coordinates,quantity,slice_time=slice_time)
    if (path=='processes_number'):
        out=load_wrf_morr_num_proc(filename,add_coordinates,quantity)
    if path=='hydrometeor':
        out=calculate_wrf_morr_path_hydrometeors(filename,slice_time=slice_time)
    if path=='phase':
        out=calculate_wrf_morr_path_phases(filename,slice_time=slice_time)
    if (path=='vaporliquid'):
        out=calculate_wrf_morr_path_vaporliquid(filename,slice_time=slice_time)
    if path=='vaporfrozen':
        out=calculate_wrf_morr_path_vaporfrozen(filename,slice_time=slice_time)
    if path=='liquidfrozen':
        out=calculate_wrf_morr_path_liquidfrozen(filename,slice_time=slice_time)
    if path=='VAPORCLOUD':
        out=calculate_wrf_morr_path_VAPORCLOUD(filename,slice_time=slice_time)
    if path=='VAPORRAIN':
        out=calculate_wrf_morr_path_VAPORRAIN(filename,slice_time=slice_time)
    if path=='VAPORICE':
        out=calculate_wrf_morr_path_VAPORICE(filename,slice_time=slice_time)   
    if path=='VAPORSNOW':
        out=calculate_wrf_morr_path_VAPORSNOW(filename,slice_time=slice_time)   
    if path=='VAPORGRAUP':
        out=calculate_wrf_morr_path_VAPORGRAUP(filename,slice_time=slice_time)    
    if path=='CLOUDRAIN':
        out=calculate_wrf_morr_path_CLOUDRAIN(filename,slice_time=slice_time)   
    if path=='CLOUDICE':
        out=calculate_wrf_morr_path_CLOUDICE(filename,slice_time=slice_time)    
    if path=='CLOUDSNOW':
        out=calculate_wrf_morr_path_CLOUDSNOW(filename,slice_time=slice_time)   
    if path=='CLOUDGRAUP':
        out=calculate_wrf_morr_path_CLOUDGRAUP(filename,slice_time=slice_time)   
    if path=='RAINICE':
        out=calculate_wrf_morr_path_RAINICE(filename,slice_time=slice_time)  
    if path=='RAINSNOW':
        out=calculate_wrf_morr_path_RAINSNOW(filename,slice_time=slice_time)
    if path=='RAINGRAUP':
        out=calculate_wrf_morr_path_RAINGRAUP(filename,slice_time=slice_time)
    if path=='ICESNOW':
        out=calculate_wrf_morr_path_ICESNOW(filename,slice_time=slice_time)
    if path=='ICEGRAUP':
        out=calculate_wrf_morr_path_ICEGRAUP(filename,slice_time=slice_time)
    if path=='SNOWGRAUP':
        out=calculate_wrf_morr_path_SNOWGRAUP(filename,slice_time=slice_time)
    else:
        print('path string unknown')
    return out   
    
def calculate_wrf_morr_latentheating(filename):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    LHREVP=loadwrfcube(filename,'LHREVP')
    LHRFRZ=loadwrfcube(filename,'LHRFRZ')
    LHRSUB=loadwrfcube(filename,'LHRSUB')
    LHR=-1*LHREVP+LHRFRZ+LHRSUB
    LHR.rename('latent heating rate')
    return LHR



def calculate_wrf_morr_path_VAPORCLOUD(filename,slice_time=slice(None)):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    #print('calculate process rates VAPOR/CLOUD')
    PCC= loadwrfcube(filename, 'PCC3D',slice_time=slice_time) 
    #PCCN=loadwrfcube(filename, 'PCCN3D')
    P_VAPORCLOUD = PCC#+PCCN
    P_VAPORCLOUD.rename('P_VAPORCLOUD')
    return P_VAPORCLOUD
    
def calculate_wrf_morr_path_VAPORRAIN(filename,slice_time=slice(None)):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    #print('calculate process rates VAPOR/RAIN')
    P_VAPORRAIN = loadwrfcube(filename, 'PRE3D',slice_time=slice_time)  #EVAP OF RAIN
    P_VAPORRAIN.rename('P_VAPORRAIN')

    return P_VAPORRAIN
    
def calculate_wrf_morr_path_VAPORICE(filename,slice_time=slice(None)):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    #print('calculate process rates VAPOR/ICE')
    PRD=loadwrfcube(filename, 'PRD3D',slice_time=slice_time)     # DEP CLOUD ICE
    EPRD=loadwrfcube(filename, 'EPRD3D',slice_time=slice_time)     # SUBLIMATION CLOUD ICE
    MNUCCD=loadwrfcube(filename, 'MNUCCD3D',slice_time=slice_time) # CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
    P_VAPORICE = PRD + EPRD + MNUCCD
    P_VAPORICE.rename('P_VAPORICE')
    return P_VAPORICE
    
def calculate_wrf_morr_path_VAPORSNOW(filename,slice_time=slice(None)):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    #print('calculate process rates VAPOR/SNOW')
    EVPMS=loadwrfcube(filename, 'EVPMS3D',slice_time=slice_time) # CHANGE Q MELTING SNOW EVAPORATING
    EPRDS=loadwrfcube(filename, 'EPRDS3D',slice_time=slice_time) #    SUBLIMATION SNOW
    P_VAPORSNOW=EVPMS+EPRDS
    P_VAPORSNOW.rename('P_VAPORSNOW')
    return P_VAPORSNOW
    
def calculate_wrf_morr_path_VAPORGRAUP(filename,slice_time=slice(None)):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    #print('calculate process rates VAPOR/GRAUPEL')
    EVPMG= loadwrfcube(filename, 'EVPMG3D',slice_time=slice_time)  # CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
    PRDG=loadwrfcube(filename,'PRDG3D',slice_time=slice_time)  # DEP OF GRAUPEL
    EPRDG=loadwrfcube(filename, 'EPRDG3D',slice_time=slice_time)  #  SUB OF GRAUPEL
    P_VAPORGRAUPEL=EVPMG+PRDG+EPRDG
    P_VAPORGRAUPEL.rename('P_VAPORGRAUPEL')
    return P_VAPORGRAUPEL
    
def calculate_wrf_morr_path_CLOUDRAIN(filename,slice_time=slice(None)):
    #Load and add up all process rates between cloud droplets and rain
    from wrfload import loadwrfcube
    #print('calculate process rates CLOUD/RAIN')
    PRA= loadwrfcube(filename, 'PRA3D',slice_time=slice_time)      # ACCRETION DROPLETS BY RAIN
    PRC= loadwrfcube(filename, 'PRC3D',slice_time=slice_time)    # AUTOCONVERSION DROPLETS  
    P_CLOUDRAIN=PRA+PRC
    P_CLOUDRAIN.rename('P_CLOUDRAIN')
    return P_CLOUDRAIN
    
    
def calculate_wrf_morr_path_CLOUDICE(filename,slice_time=slice(None)):
    #Load and add up all process rates between ckoud droplets and cloud ice
    from wrfload import loadwrfcube
    # print('calculate process rates CLOUD/ICE')
    PSACWI=loadwrfcube(filename, 'PSACWI3D',slice_time=slice_time)   # CHANGE Q DROPLET ACCRETION BY CLOUD ICE
    QMULTS=loadwrfcube(filename, 'QMULTS3D',slice_time=slice_time)  # CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
    QMULTG=loadwrfcube(filename, 'QMULTG3D',slice_time=slice_time)   # CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
    P_CLOUDICE =PSACWI+QMULTS+QMULTG
    P_CLOUDICE.rename('P_CLOUDICE')
    return P_CLOUDICE
    
    
def calculate_wrf_morr_path_CLOUDSNOW(filename,slice_time=slice(None)):
    #Load and add up all process rates between cloud droplets and snow
    from wrfload import loadwrfcube
    # print('calculate process rates CLOUD/ICE')
    PSACWS= loadwrfcube(filename, 'PSACWS3D',slice_time=slice_time)      # CHANGE Q DROPLET ACCRETION BY SNOW
    P_CLOUDSNOW=PSACWS
    P_CLOUDSNOW.rename('P_CLOUDSNOW')
    return P_CLOUDSNOW

def calculate_wrf_morr_path_CLOUDGRAUP(filename,slice_time=slice(None)):
    #Load and add up all process rates between cloud droplets and graupel
    from wrfload import loadwrfcube
    # print('calculate process rates CLOUD/GRAUP')
    PSACWG=loadwrfcube(filename, 'PSACWG3D',slice_time=slice_time)      #  CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
    PGSACW=loadwrfcube(filename, 'PGSACW3D',slice_time=slice_time)   # CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
    P_CLOUDGRAUP = PGSACW + PSACWG
    P_CLOUDGRAUP.rename('P_CLOUDGRAUP')
    return P_CLOUDGRAUP


def calculate_wrf_morr_path_RAINICE(filename,slice_time=slice(None)):
    #Load and add up all process rates between rain and cloud ice
    from wrfload import loadwrfcube
    # print('calculate process rates RAIN/ICE')
    QMULTR=loadwrfcube(filename, 'QMULTR3D',slice_time=slice_time)      # CHANGE Q DUE TO ICE RAIN/SNOW
    QMULTRG=loadwrfcube(filename, 'QMULTRG3D',slice_time=slice_time)                         # CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL
    P_RAINICE = QMULTR+QMULTRG
    P_RAINICE.rename('P_RAINICE')
    return P_RAINICE    

def calculate_wrf_morr_path_RAINSNOW(filename,slice_time=slice(None)):
    #Load and add up all process rates between rain and snow
    from wrfload import loadwrfcube
    # print('calculate process rates RAIN/SNOW')
    PSMLT=loadwrfcube(filename, 'PSMLT3D',slice_time=slice_time)     # CHANGE Q MELTING SNOW TO RAIN
    PIACRS=loadwrfcube(filename, 'PIACRS3D',slice_time=slice_time)                           #CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
    P_RAINSNOW =PSMLT+PIACRS
    P_RAINSNOW.rename('P_RAINSNOW')
    return P_RAINSNOW   

def calculate_wrf_morr_path_RAINGRAUP(filename,slice_time=slice(None)):
    #Load and add up all process rates between rain and graupel
    from wrfload import loadwrfcube
    # print('calculate process rates RAIN/GRAUPEL')
    MNUCCR=loadwrfcube(filename, 'MNUCCR3D',slice_time=slice_time)     # CHANGE Q DUE TO CONTACT FREEZ RAIN
    PIACR=loadwrfcube(filename, 'PIACR3D',slice_time=slice_time)      # CHANGE QR, ICE-RAIN COLLECTION
    PRACG=loadwrfcube(filename, 'PRACG3D',slice_time=slice_time)      #CHANGE IN Q COLLECTION RAIN BY GRAUPEL
    PGRACS=loadwrfcube(filename, 'PGRACS3D',slice_time=slice_time) # CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
    PGMLT=loadwrfcube(filename, 'PGMLT3D',slice_time=slice_time) #  CHANGE Q MELTING OF GRAUPEL
    P_RAINGRAUP =MNUCCR+PIACR+PRACG+PGRACS+PGMLT
    P_RAINGRAUP.rename('P_RAINGRAUP')
    return P_RAINGRAUP

def calculate_wrf_morr_path_ICESNOW(filename,slice_time=slice(None)):
    #Load and add up all process rates between cloud ice and snow
    from wrfload import loadwrfcube
    # print('calculate process rates ICE/SNOW')
    PRAI = loadwrfcube(filename, 'PRAI3D',slice_time=slice_time)      # CHANGE Q ACCRETION CLOUD ICE BY SNOW
    PRCI=loadwrfcube(filename, 'PRCI3D',slice_time=slice_time)      # CHANGE Q AUTOCONVERSIN CLOUD ICE TO SNOW
    PRACIS=loadwrfcube(filename, 'PRACIS3D',slice_time=slice_time)     # CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
    P_ICESNOW = PRAI + PRCI + PRACIS
    P_ICESNOW.rename('P_ICESNOW')
    return P_ICESNOW  
    
def calculate_wrf_morr_path_ICEGRAUP(filename,slice_time=slice(None)):
    #Load and add up all process rates between cloud ice and graupel
    from wrfload import loadwrfcube
    # print('calculate process rates ICE/GRAUPEL')
    PRACI=loadwrfcube(filename, 'PRACI3D',slice_time=slice_time)     # CHANGE QI, ICE-RAIN COLLECTION
    P_ICEGRAUP = PRACI
    P_ICEGRAUP .rename('P_ICEGRAUP ')
    return P_ICEGRAUP
    
def calculate_wrf_morr_path_SNOWGRAUP(filename,slice_time=slice(None)):
    #Load and add up all process rates between snow and graupel
    from wrfload import loadwrfcube
    # print('calculate process rates SNOW/GRAUPEL')
    P_SNOWGRAUP = 0*loadwrfcube(filename, 'PRACI3D',slice_time=slice_time)   # Dummy zeros, since no pathway process found yet
    P_SNOWGRAUP.rename('P_SNOWGRAUP')
    return P_SNOWGRAUP  
    
        
def calculate_wrf_morr_path_vaporliquid(filename,slice_time=slice(None)):
    #Load and add up all process rates between ice phase and water vapour:
    #print('calculate processes deposition/sublimation')
    PVAPORCLOUD=calculate_wrf_morr_path_VAPORCLOUD(filename,slice_time=slice_time)
    PVAPORRAIN=calculate_wrf_morr_path_VAPORRAIN(filename,slice_time=slice_time)  
    P_vaporliquid=  PVAPORCLOUD + PVAPORRAIN
    P_vaporliquid.rename('PVAPORLIQUID')
    return P_vaporliquid
    
def calculate_wrf_morr_path_vaporfrozen(filename,slice_time=slice(None)):
    #Load and add up all process rates between ice phase and water vapour:
    #print('calculate processes deposition/sublimation')
    PVAPORICE=calculate_wrf_morr_path_VAPORICE(filename,slice_time=slice_time)
    PVAPORSNOW=calculate_wrf_morr_path_VAPORSNOW(filename,slice_time=slice_time)
    PVAPORGRAUP=calculate_wrf_morr_path_VAPORGRAUP(filename,slice_time=slice_time)   
    P_vaporfrozen=PVAPORICE+PVAPORSNOW+PVAPORGRAUP
    P_vaporfrozen.rename('PVAPORFROZEN')
    return P_vaporfrozen    
    
def calculate_wrf_morr_path_liquidfrozen(filename,slice_time=slice(None)):
    #Load and add up all process rates between frozen and liquid phase
    #print('calculate processes freezing/melting')
    PCLOUDICE=calculate_wrf_morr_path_CLOUDICE(filename,slice_time=slice_time)
    PRAINICE=calculate_wrf_morr_path_RAINICE(filename,slice_time=slice_time)
    PCLOUDSNOW=calculate_wrf_morr_path_CLOUDSNOW(filename,slice_time=slice_time)
    PRAINSNOW=calculate_wrf_morr_path_RAINSNOW(filename,slice_time=slice_time)
    PRAINGRAUP=calculate_wrf_morr_path_RAINGRAUP(filename,slice_time=slice_time)
    PCLOUDGRAUP=calculate_wrf_morr_path_CLOUDGRAUP(filename,slice_time=slice_time)
    P_liquidfrozen=PCLOUDICE+PRAINICE+PCLOUDSNOW+PRAINSNOW+PRAINGRAUP+PCLOUDGRAUP
    P_liquidfrozen.rename('PLIQUIDFROZEN')

    return P_liquidfrozen    

def sum_cubes(filename,name,list_names):
    from wrfload import loadwrfcube
    P_out=loadwrfcube(filename, list_names[0])
    for name_i in listnames[1:]:
        P_out=P_out+loadwrfcube(filename, name_i)
    P_out.rename(name)    
    return P_out

    
    

def calculate_wrf_thom_path_VAPORCLOUD(filename):
    list_names=['']
    name='P_VAPORCLOUD'
    return sum_cubes(filename,name,list_names)
    
def calculate_wrf_thom_path_VAPORRAIN(filename):
    #Load and add up all process rates between water vapour and cloud droplets
    list_names=['']
    name='P_VAPORRAIN'
    return sum_cubes(filename,name,list_names)
    
def calculate_wrf_thom_path_VAPORICE(filename):
    #Load and add up all process rates between water vapour and cloud droplets
    list_names=['']
    name='P_VAPORICE'
    return sum_cubes(filename,name,list_names)
    
def calculate_wrf_thom_path_VAPORSNOW(filename):
    list_names=['']
    name='P_VAPORSNOW'
    return sum_cubes(filename,name,list_names)
    
def calculate_wrf_thom_path_VAPORGRAUP(filename):
    #Load and add up all process rates between water vapour and cloud droplets
    list_names=['']
    name='P_VAPORGRAUPEL'
    return sum_cubes(filename,name,list_names)
    
def calculate_wrf_thom_path_CLOUDRAIN(filename):
     #Load and add up all process rates between water vapour and cloud droplets
    list_names=['']
    name='P_CLOUDRAIN'
    return sum_cubes(filename,name,list_names)
    
    
def calculate_wrf_thom_path_CLOUDICE(filename):
    #Load and add up all process rates between ckoud droplets and cloud ice
    list_names=['']
    name='P_CLOUDICE'
    return sum_cubes(filename,name,list_names)
    
    
    
def calculate_wrf_thom_path_CLOUDSNOW(filename):
    #Load and add up all process rates between cloud droplets and snow
    list_names=['']
    name='P_CLOUDSNOW'
    return sum_cubes(filename,name,list_names)

def calculate_wrf_thom_path_CLOUDGRAUP(filename):
    #Load and add up all process rates between cloud droplets and graupel
    list_names=['']
    name='P_CLOUDGRAUP'
    return sum_cubes(filename,name,list_names)


def calculate_wrf_thom_path_RAINICE(filename):
    #Load and add up all process rates between rain and cloud ice
    list_names=['']
    name='P_RAINICE'
    return sum_cubes(filename,name,list_names)
def calculate_wrf_thom_path_RAINSNOW(filename):
    #Load and add up all process rates between rain and snow
    list_names=['']
    name='P_RAINICE'
    return sum_cubes(filename,name,list_names)

def calculate_wrf_thom_path_RAINGRAUP(filename):
    #Load and add up all process rates between rain and graupel
    list_names=['']
    name='P_RAINGRAUP'
    return sum_cubes(filename,name,list_names)

def calculate_wrf_thom_path_ICESNOW(filename):
    #Load and add up all process rates between cloud ice and snow
    list_names=['']
    name='P_ICESNOW'
    return sum_cubes(filename,name,list_names)
    
def calculate_wrf_thom_path_ICEGRAUP(filename):
    #Load and add up all process rates between cloud ice and graupel
    list_names=['']
    name='P_ICEGRAUP'
    return sum_cubes(filename,name,list_names)
    
def calculate_wrf_thom_path_SNOWGRAUP(filename):
    #Load and add up all process rates between snow and graupel
    list_names=['']
    name='P_SNOWGRAUP'
    return sum_cubes(filename,name,list_names)
    
    
def calculate_wrf_thom_path_vaporliquid(filename):
    #Load and add up all process rates between ice phase and water vapour:
    #print('calculate processes deposition/sublimation')
    PVAPORCLOUD=calculate_wrf_thom_path_VAPORCLOUD(filename)
    PVAPORRAIN=calculate_wrf_thom_path_VAPORRAIN(filename)  
    P_vaporliquid=  PVAPORCLOUD + PVAPORRAIN
    P_vaporliquid.rename('PVAPORLIQUID')
    return sum_cubes(filename,name,list_names)
    
def calculate_wrf_thom_path_vaporfrozen(filename):
    #Load and add up all process rates between ice phase and water vapour:
    #print('calculate processes deposition/sublimation')
    PVAPORICE=calculate_wrf_thom_path_VAPORICE(filename)
    PVAPORSNOW=calculate_wrf_thom_path_VAPORSNOW(filename)
    PVAPORGRAUP=calculate_wrf_thom_path_VAPORGRAUP(filename)   
    P_vaporfrozen=PVAPORICE+PVAPORSNOW+PVAPORGRAUP
    P_vaporfrozen.rename('PVAPORFROZEN')
    return P_vaporfrozen    
    
def calculate_wrf_thom_path_liquidfrozen(filename):
    #Load and add up all process rates between frozen and liquid phase
    #print('calculate processes freezing/melting')
    PCLOUDICE=calculate_wrf_thom_path_CLOUDICE(filename)
    PRAINICE=calculate_wrf_thom_path_RAINICE(filename)
    PCLOUDSNOW=calculate_wrf_thom_path_CLOUDSNOW(filename)
    PRAINSNOW=calculate_wrf_thom_path_RAINSNOW(filename)
    PRAINGRAUP=calculate_wrf_thom_path_RAINGRAUP(filename)
    PCLOUDGRAUP=calculate_wrf_thom_path_CLOUDGRAUP(filename)
    P_liquidfrozen=PCLOUDICE+PRAINICE+PCLOUDSNOW+PRAINSNOW+PRAINGRAUP+PCLOUDGRAUP
    P_liquidfrozen.rename('PLIQUIDFROZEN')

    return P_liquidfrozen    


    
    
    
