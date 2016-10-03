# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 12:30:47 2016


Functions for the analysis of WRF microphysical pathways
@author: heikenfeld
"""

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
Proclist_Morr_mass_load=['PRD',
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

def load_wrf_morr_mass_proc(filename):
    from wrfload import loadwrfcube
    import numpy as np
    Dict={}
    for process in Proclist_Morr_mass_load:
        print(process)
        if process=='PCC':
            cube=loadwrfcube(filename,process+'3D')
            cube1=cube[:]
            cube1.data=np.clip(cube1.data,a_min=-np.inf,a_max=0)            
            cube1.rename('PCC')
            Dict['PCC']=cube1
            cube2=cube[:]
            cube2.data=np.clip(cube1.data,a_min=0,a_max=np.inf)
            cube2.rename('EPCC')
            Dict['EPCC']=cube2

        else:
            cube=loadwrfcube(filename,process+'3D')
            cube.rename('process')
            #Cubelist.append(cube)
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
    
def load_wrf_morr_num_proc(filename):
    from wrfload import loadwrfcube

 
    #Cubelist=[]
    Dict={}
    for process in Proclist_Morr_number:
        if process=='NSMLTR':
                cube=loadwrfcube(filename,process)
        else:
            cube=loadwrfcube(filename,process+'3D')
            cube.rename('process')
        #Cubelist.append(cube)
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
    
def calculate_wrf_morr_path_hydrometeors(filename):
    Dict={}
    #Cubelist=[]
    for path in Hydropath_list:
        cube=calculate_wrf_morr_path(filename,path)
        #Cubelist.append(cube)
        Dict[path]=cube
    #return Cubelist
    return Dict
    
Phasepath_list=['vaporliquid','vaporfrozen','liquidfrozen']

def calculate_wrf_morr_path_phases(filename):
    Dict={}
    #Cubelist=[]
    for path in Phasepath_list:
        cube=calculate_wrf_morr_path(filename,path)
        #Cubelist.append(cube)
        Dict[path]=cube
    #return Cubelist
    return Dict


        
def calculate_wrf_morr_path(filename,path):
    if (path=='processes_mass'):
        out=load_wrf_morr_mass_proc(filename)
    if (path=='processes_number'):
        out=load_wrf_morr_num_proc(filename)
    if path=='hydrometeor':
        out=calculate_wrf_morr_path_hydrometeors(filename)
    if path=='phase':
        out=calculate_wrf_morr_path_phases(filename)
    if (path=='vaporliquid'):
        out=calculate_wrf_morr_path_vaporliquid(filename)
    if path=='vaporfrozen':
        out=calculate_wrf_morr_path_vaporfrozen(filename)
    if path=='liquidfrozen':
        out=calculate_wrf_morr_path_liquidfrozen(filename)
    if path=='VAPORCLOUD':
        out=calculate_wrf_morr_path_VAPORCLOUD(filename)
    if path=='VAPORRAIN':
        out=calculate_wrf_morr_path_VAPORRAIN(filename)
    if path=='VAPORICE':
        out=calculate_wrf_morr_path_VAPORICE(filename)   
    if path=='VAPORSNOW':
        out=calculate_wrf_morr_path_VAPORSNOW(filename)   
    if path=='VAPORGRAUP':
        out=calculate_wrf_morr_path_VAPORGRAUP(filename)    
    if path=='CLOUDRAIN':
        out=calculate_wrf_morr_path_CLOUDRAIN(filename)   
    if path=='CLOUDICE':
        out=calculate_wrf_morr_path_CLOUDICE(filename)    
    if path=='CLOUDSNOW':
        out=calculate_wrf_morr_path_CLOUDSNOW(filename)   
    if path=='CLOUDGRAUP':
        out=calculate_wrf_morr_path_CLOUDGRAUP(filename)   
    if path=='RAINICE':
        out=calculate_wrf_morr_path_RAINICE(filename)  
    if path=='RAINSNOW':
        out=calculate_wrf_morr_path_RAINSNOW(filename)
    if path=='RAINGRAUP':
        out=calculate_wrf_morr_path_RAINGRAUP(filename)
    if path=='ICESNOW':
        out=calculate_wrf_morr_path_ICESNOW(filename)
    if path=='ICEGRAUP':
        out=calculate_wrf_morr_path_ICEGRAUP(filename)
    if path=='SNOWGRAUP':
        out=calculate_wrf_morr_path_SNOWGRAUP(filename)
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



def calculate_wrf_morr_path_VAPORCLOUD(filename):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    #print('calculate process rates VAPOR/CLOUD')
    PCC= loadwrfcube(filename, 'PCC3D') 
    #PCCN=loadwrfcube(filename, 'PCCN3D')
    P_VAPORCLOUD = PCC#+PCCN
    P_VAPORCLOUD.rename('P_VAPORCLOUD')
    return P_VAPORCLOUD
    
def calculate_wrf_morr_path_VAPORRAIN(filename):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    #print('calculate process rates VAPOR/RAIN')
    P_VAPORRAIN = loadwrfcube(filename, 'PRE3D')  #EVAP OF RAIN
    P_VAPORRAIN.rename('P_VAPORRAIN')

    return P_VAPORRAIN
    
def calculate_wrf_morr_path_VAPORICE(filename):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    #print('calculate process rates VAPOR/ICE')
    PRD=loadwrfcube(filename, 'PRD3D')     # DEP CLOUD ICE
    EPRD=loadwrfcube(filename, 'EPRD3D')     # SUBLIMATION CLOUD ICE
    MNUCCD=loadwrfcube(filename, 'MNUCCD3D') # CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
    P_VAPORICE = PRD + EPRD + MNUCCD
    P_VAPORICE.rename('P_VAPORICE')
    return P_VAPORICE
    
def calculate_wrf_morr_path_VAPORSNOW(filename):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    #print('calculate process rates VAPOR/SNOW')
    EVPMS=loadwrfcube(filename, 'EVPMS3D') # CHANGE Q MELTING SNOW EVAPORATING
    EPRDS=loadwrfcube(filename, 'EPRDS3D') #    SUBLIMATION SNOW
    P_VAPORSNOW=EVPMS+EPRDS
    P_VAPORSNOW.rename('P_VAPORSNOW')
    return P_VAPORSNOW
    
def calculate_wrf_morr_path_VAPORGRAUP(filename):
    #Load and add up all process rates between water vapour and cloud droplets
    from wrfload import loadwrfcube
    #print('calculate process rates VAPOR/GRAUPEL')
    EVPMG= loadwrfcube(filename, 'EVPMG3D')  # CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
    PRDG=loadwrfcube(filename,'PRDG3D')  # DEP OF GRAUPEL
    EPRDG=loadwrfcube(filename, 'EPRDG3D')  #  SUB OF GRAUPEL
    P_VAPORGRAUPEL=EVPMG+PRDG+EPRDG
    P_VAPORGRAUPEL.rename('P_VAPORGRAUPEL')
    return P_VAPORGRAUPEL
    
def calculate_wrf_morr_path_CLOUDRAIN(filename):
    #Load and add up all process rates between cloud droplets and rain
    from wrfload import loadwrfcube
    #print('calculate process rates CLOUD/RAIN')
    PRA= loadwrfcube(filename, 'PRA3D')      # ACCRETION DROPLETS BY RAIN
    PRC= loadwrfcube(filename, 'PRC3D')    # AUTOCONVERSION DROPLETS  
    P_CLOUDRAIN=PRA+PRC
    P_CLOUDRAIN.rename('P_CLOUDRAIN')
    return P_CLOUDRAIN
    
    
def calculate_wrf_morr_path_CLOUDICE(filename):
    #Load and add up all process rates between ckoud droplets and cloud ice
    from wrfload import loadwrfcube
    # print('calculate process rates CLOUD/ICE')
    PSACWI=loadwrfcube(filename, 'PSACWI3D')   # CHANGE Q DROPLET ACCRETION BY CLOUD ICE
    QMULTS=loadwrfcube(filename, 'QMULTS3D')  # CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
    QMULTG=loadwrfcube(filename, 'QMULTG3D')   # CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
    P_CLOUDICE =PSACWI+QMULTS+QMULTG
    P_CLOUDICE.rename('P_CLOUDICE')
    return P_CLOUDICE
    
    
def calculate_wrf_morr_path_CLOUDSNOW(filename):
    #Load and add up all process rates between cloud droplets and snow
    from wrfload import loadwrfcube
    # print('calculate process rates CLOUD/ICE')
    PSACWS= loadwrfcube(filename, 'PSACWS3D')      # CHANGE Q DROPLET ACCRETION BY SNOW
    P_CLOUDSNOW=PSACWS
    P_CLOUDSNOW.rename('P_CLOUDSNOW')
    return P_CLOUDSNOW

def calculate_wrf_morr_path_CLOUDGRAUP(filename):
    #Load and add up all process rates between cloud droplets and graupel
    from wrfload import loadwrfcube
    # print('calculate process rates CLOUD/GRAUP')
    PSACWG=loadwrfcube(filename, 'PSACWG3D')      #  CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
    PGSACW=loadwrfcube(filename, 'PGSACW3D')   # CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
    P_CLOUDGRAUP = PGSACW + PSACWG
    P_CLOUDGRAUP.rename('P_CLOUDGRAUP')
    return P_CLOUDGRAUP


def calculate_wrf_morr_path_RAINICE(filename):
    #Load and add up all process rates between rain and cloud ice
    from wrfload import loadwrfcube
    # print('calculate process rates RAIN/ICE')
    QMULTR=loadwrfcube(filename, 'QMULTR3D')      # CHANGE Q DUE TO ICE RAIN/SNOW
    QMULTRG=loadwrfcube(filename, 'QMULTRG3D')                         # CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL
    P_RAINICE = QMULTR+QMULTRG
    P_RAINICE.rename('P_RAINICE')
    return P_RAINICE    

def calculate_wrf_morr_path_RAINSNOW(filename):
    #Load and add up all process rates between rain and snow
    from wrfload import loadwrfcube
    # print('calculate process rates RAIN/SNOW')
    PSMLT=loadwrfcube(filename, 'PSMLT3D')     # CHANGE Q MELTING SNOW TO RAIN
    PIACRS=loadwrfcube(filename, 'PIACRS3D')                           #CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
    P_RAINSNOW =PSMLT+PIACRS
    P_RAINSNOW.rename('P_RAINSNOW')
    return P_RAINSNOW   

def calculate_wrf_morr_path_RAINGRAUP(filename):
    #Load and add up all process rates between rain and graupel
    from wrfload import loadwrfcube
    # print('calculate process rates RAIN/GRAUPEL')
    MNUCCR=loadwrfcube(filename, 'MNUCCR3D')     # CHANGE Q DUE TO CONTACT FREEZ RAIN
    PIACR=loadwrfcube(filename, 'PIACR3D')      # CHANGE QR, ICE-RAIN COLLECTION
    PRACG=loadwrfcube(filename, 'PRACG3D')      #CHANGE IN Q COLLECTION RAIN BY GRAUPEL
    PGRACS=loadwrfcube(filename, 'PGRACS3D') # CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
    PGMLT=loadwrfcube(filename, 'PGMLT3D') #  CHANGE Q MELTING OF GRAUPEL
    P_RAINGRAUP =MNUCCR+PIACR+PRACG+PGRACS+PGMLT
    P_RAINGRAUP.rename('P_RAINGRAUP')
    return P_RAINGRAUP

def calculate_wrf_morr_path_ICESNOW(filename):
    #Load and add up all process rates between cloud ice and snow
    from wrfload import loadwrfcube
    # print('calculate process rates ICE/SNOW')
    PRAI = loadwrfcube(filename, 'PRAI3D')      # CHANGE Q ACCRETION CLOUD ICE BY SNOW
    PRCI=loadwrfcube(filename, 'PRCI3D')      # CHANGE Q AUTOCONVERSIN CLOUD ICE TO SNOW
    PRACIS=loadwrfcube(filename, 'PRACIS3D')     # CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
    P_ICESNOW = PRAI + PRCI + PRACIS
    P_ICESNOW.rename('P_ICESNOW')
    return P_ICESNOW  
    
def calculate_wrf_morr_path_ICEGRAUP(filename):
    #Load and add up all process rates between cloud ice and graupel
    from wrfload import loadwrfcube
    # print('calculate process rates ICE/GRAUPEL')
    PRACI=loadwrfcube(filename, 'PRACI3D')     # CHANGE QI, ICE-RAIN COLLECTION
    P_ICEGRAUP = PRACI
    P_ICEGRAUP .rename('P_ICEGRAUP ')
    return P_ICEGRAUP
    
def calculate_wrf_morr_path_SNOWGRAUP(filename):
    #Load and add up all process rates between snow and graupel
    from wrfload import loadwrfcube
    # print('calculate process rates SNOW/GRAUPEL')
    P_SNOWGRAUP = 0*loadwrfcube(filename, 'PRACI3D')   # Dummy zeros, since no pathway process found yet
    P_SNOWGRAUP.rename('P_SNOWGRAUP')
    return P_SNOWGRAUP  
    
        
def calculate_wrf_morr_path_vaporliquid(filename):
    #Load and add up all process rates between ice phase and water vapour:
    #print('calculate processes deposition/sublimation')
    PVAPORCLOUD=calculate_wrf_morr_path_VAPORCLOUD(filename)
    PVAPORRAIN=calculate_wrf_morr_path_VAPORRAIN(filename)  
    P_vaporliquid=  PVAPORCLOUD + PVAPORRAIN
    P_vaporliquid.rename('PVAPORLIQUID')
    return P_vaporliquid
    
def calculate_wrf_morr_path_vaporfrozen(filename):
    #Load and add up all process rates between ice phase and water vapour:
    #print('calculate processes deposition/sublimation')
    PVAPORICE=calculate_wrf_morr_path_VAPORICE(filename)
    PVAPORSNOW=calculate_wrf_morr_path_VAPORSNOW(filename)
    PVAPORGRAUP=calculate_wrf_morr_path_VAPORGRAUP(filename)   
    P_vaporfrozen=PVAPORICE+PVAPORSNOW+PVAPORGRAUP
    P_vaporfrozen.rename('PVAPORFROZEN')
    return P_vaporfrozen    
    
def calculate_wrf_morr_path_liquidfrozen(filename):
    #Load and add up all process rates between frozen and liquid phase
    #print('calculate processes freezing/melting')
    PCLOUDICE=calculate_wrf_morr_path_CLOUDICE(filename)
    PRAINICE=calculate_wrf_morr_path_RAINICE(filename)
    PCLOUDSNOW=calculate_wrf_morr_path_CLOUDSNOW(filename)
    PRAINSNOW=calculate_wrf_morr_path_RAINSNOW(filename)
    PRAINGRAUP=calculate_wrf_morr_path_RAINGRAUP(filename)
    PCLOUDGRAUP=calculate_wrf_morr_path_CLOUDGRAUP(filename)
    P_liquidfrozen=PCLOUDICE+PRAINICE+PCLOUDSNOW+PRAINSNOW+PRAINGRAUP+PCLOUDGRAUP
    P_liquidfrozen.rename('PLIQUIDFROZEN')

    return P_liquidfrozen    
    
    
