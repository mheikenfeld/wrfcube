from collections import defaultdict


def split_sign_variable(filename,variable,name_neg,name_pos,add_coordinates=None,constraint=None):
   from wrfload import loadwrfcube
   cube=loadwrfcube(filename,variable,add_coordinates=add_coordinates,constraint=constraint)
   #dict_out={}
   list_out=[]
   #dict_out[name_neg]=get_process_neg(cube,cube,name_neg)
   #dict_out[name_pos]=get_process_pos(cube,cube,name_pos)
   list_out.append(get_variable_neg(cube,name_neg))
   list_out.append( get_variable_pos(cube,name_pos))
   #return dict_out
   return list_out

def get_variable_pos(cube,name_neg):
   import numpy as np
   cube_neg=cube[:]
   cube_neg.data=np.clip(cube.data,a_max=np.inf,a_min=0)
   cube_neg.rename(name_neg)
   return cube_neg

def get_variable_neg(cube,name_pos):
   import numpy as np
   cube_pos=cube[:]
   cube_pos.data=np.abs(np.clip(cube.data,a_min=-np.inf,a_max=0))
   cube_pos.rename(name_pos)

   return cube_pos
   


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
         'PRG_IHM',   #  Graupel->Ice
         'PRW_IMI',   #  Ice -> Water
         'PRI_WFI'    #  Water -> Ice
          ]
          
          
thompson_processes_mass_split=defaultdict(dict)
thompson_processes_mass_split['PRW_VCD']=['PRW_VCD','E_PRW_VCD']
thompson_processes_mass_split['PRR_RCS']=['PRR_RCS','E_PRR_RCS']
thompson_processes_mass_split['PRR_RCG']=['PRR_RCG','E_PRR_RCG']
thompson_processes_mass_split['PRI_IDE']=['PRI_IDE','E_PRI_IDE']
thompson_processes_mass_split['PRS_RCS']=['PRS_RCS','E_PRS_RCS']
thompson_processes_mass_split['PRS_SDE']=['PRS_SDE','E_PRS_SDE']
thompson_processes_mass_split['PRG_GDE']=['PRG_GDE','E_PRG_GDE']
thompson_processes_mass_split['PRG_RCG']=['PRG_RCG','E_PRG_RCG']

List_Processes_Thompson_Mass_signed=list(List_Processes_Thompson_Mass).extend(['E_PRW_VCD,E_PRR_RCS','E_PRR_RCG','E_PRI_IDE','E_PRS_RCS','E_PRS_SDE','E_PRG_GDE','E_PRG_RCG'])


#Processes_Thompson_Mass_colors={}
#Processes_Thompson_Mass_colors['PRW_VCD']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRV_REV']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRR_WAU']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRR_RCW']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRR_RCS']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRR_RCG']='slategrey'   #  Rain->Graupel
#Processes_Thompson_Mass_colors['PRR_GML']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRR_RCI']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRI_INU']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRI_IHM']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRI_WFZ']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRI_RFZ']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRI_IDE']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRI_RCI']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRI_IHA']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRS_IAU']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRS_SCI']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRS_RCS']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRS_SCW']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRS_SDE']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRS_IHM']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRS_IDE']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRG_SCW']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRG_RFZ']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRG_GDE']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRG_GCW']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRG_RCI']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRG_RCS']='slategrey'   #  Vapor->Water
#Processes_Thompson_Mass_colors['PRG_RCG']='slategrey'   #  Rain->Graupel
#Processes_Thompson_Mass_colors['PRG_IHM']='slategrey'   #  Graupel->Ice




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
         'PND_GCD',
         'PNW_IMI',   #  Ice -> Water
         'PNI_WFI'    #  Water -> Ice

                              ]





#def load_wrf_thom_mass_proc(filename,add_coordinates=None,quantity='volume'):
#    from wrfload import loadwrfcube, derivewrfcube
#    from iris.cube import CubeList
#    #Dict={}
#    cubelist_out=CubeList()
#
#    if quantity=='volume':
#        rho=derivewrfcube(filename,'density')
#
#    for i_process,process in enumerate(List_Processes_Thompson_Mass):
#        #print(process)
#        if (i_process==0):
#            cube=loadwrfcube(filename,process,add_coordinates=add_coordinates)
#            cube.rename(process)
#            #Cubelist.append(cube)
#            #Dict[process]=cube
#            if add_coordinates=='pz':
#               z_coord=cube.coord('geopotential')
#               p_coord=cube.coord('pressure')
#        else:
#            cube=loadwrfcube(filename,process)
#            cube.rename(process)
#            if add_coordinates=='pz':
#               cube.add_aux_coord(z_coord,(0,1,2,3))
#               cube.add_aux_coord(p_coord,(0,1,2,3))
#        if quantity=='volume':
#            cube=cube*rho
#        #Dict[process]=cube
#        cubelist_out.append(cube)
#
#    return cubelist_out
#    #return Dict
    
    
#def load_wrf_thom_mass_proc_signed(filename,add_coordinates=None,quantity='volume'):
#    from wrfload import loadwrfcube, derivewrfcube
#    from iris.cube import CubeList
#    #Dict={}
#    cubelist_out=CubeList()
#
#    if quantity=='volume':
#        rho=derivewrfcube(filename,'density')
#
#    for i_process,process in enumerate(List_Processes_Thompson_Mass):
#        #print(process)
#        if (i_process==0):
#            cube=loadwrfcube(filename,process,add_coordinates=add_coordinates)
#            cube.rename(process)
#            #Cubelist.append(cube)
#            #Dict[process]=cube
#            if add_coordinates=='pz':
#               z_coord=cube.coord('geopotential')
#               p_coord=cube.coord('pressure')
#        else:
#            cube=loadwrfcube(filename,process)
#            cube.rename(process)
#            if add_coordinates=='pz':
#               cube.add_aux_coord(z_coord,(0,1,2,3))
#               cube.add_aux_coord(p_coord,(0,1,2,3))
#        if quantity=='volume':
#            cube=cube*rho
#        #Dict[process]=cube
#        cubelist_out.append(cube)
#
#    return cubelist_out
#    #return Dict
#    
#    
    
#    
#
#def load_wrf_thom_number_proc(filename,add_coordinates=None,quantity='volume'):
#    from wrfload import loadwrfcube, derivewrfcube
#    from iris.cube import CubeList
#
#    #Dict={}
#    cubelist_out=CubeList()
#    if quantity=='volume':
#        rho=derivewrfcube(filename,'density')
#
#    for i_process,process in enumerate(List_Processes_Thompson_Number):
#        #print(process)
#        if (i_process==0):
#            cube=loadwrfcube(filename,process,add_coordinates=add_coordinates)
#            cube.rename(process)
#            #Cubelist.append(cube)
#            #Dict[process]=cube
#            if add_coordinates=='pz':
#               z_coord=cube.coord('geopotential')
#               p_coord=cube.coord('pressure')
#        else:
#            cube=loadwrfcube(filename,process)
#            cube.rename(process)
#            if add_coordinates=='pz':
#                cube.add_aux_coord(z_coord,(0,1,2,3))
#                cube.add_aux_coord(p_coord,(0,1,2,3))
#        #Cubelist.append(cube)
#        if quantity=='volume':
#            cube=cube*rho
#
#        #Dict[process]=cube
#        cubelist_out.append(cube)
#
#    #return cubelist_out
#    return Dict




#Morrison microphysics

#Morrison microphysics
morrison_processes_mass=[
    'PRD',
    'EPRD',
    'PRDG',
    'EPRDG',
    'PRDS',
    'EPRDS',
    'PRE',
    'PRA',
    'PRC',
    'PCC',
    'PCCN',
    'PSMLT',
    'EVPMS',
    'QMULTS',
    'QMULTR',
    'PRACS',
    'PSACWG',
    'PGSACW',
    'PGRACS',
    'EVPMG',
    'PGMLT',
    'PRACI',
    'PIACRS',
    'PRACIS',
    'PRACG',
    'QMULTG',
    'MNUCCR',
    'MNUCCC',
    'MNUCCD',
    'QMULTRG',
    'PRAI',
    'PRCI',
    'PSACWS',
    'PIACR',
    'PSACWI',
    'PSACR',
    'QICF', #instantaneous processes (addes by BAW)             
    'QGRF',      
    'QNIRF',
    'QIIM'
    ]
    
Proclist_Morr_mass=list(morrison_processes_mass)

morrison_processes_mass_split=defaultdict(dict)
morrison_processes_mass_split['PCC']=['EPCC','PCC']
morrison_processes_mass_split['PRCI']=['EPRCI','PRCI']

Proclist_Morr_mass_signed=list(Proclist_Morr_mass).extend(['EPCC','EPRCI'])

morrison_processes_number=[
    'NSUBC',
    'NSUBI',
    'NSUBS',
    'NSUBR',
    'NNUCCC',
    'NNUCCD',
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
    'NPSACWG',
    'NSCNG',
    'NGRACS',
    'NGMLTG',
    'NGMLTR',
    'NSUBG',
    'NMULTG',
    'NMULTRG,',
    'NIACRS', 
    'NICF', #instantaneous processes (addes by BAW)          
    'NGRF',      
    'NSRF',
    'NIIM'
]

#Morr_Processes_signed_colors=defaultdict(dict)
#Morr_Processes_signed_colors['PRD']='lightseagreen'
#Morr_Processes_signed_colors['PRE']='red'
#Morr_Processes_signed_colors['PRDS']='darkorange'
#Morr_Processes_signed_colors['PRA']='darkred'
#Morr_Processes_signed_colors['PRC']='orange'
#Morr_Processes_signed_colors['PCC']='magenta'
#Morr_Processes_signed_colors['PCCN']='blue'
#Morr_Processes_signed_colors['PSMLT']='darkblue'
#Morr_Processes_signed_colors['EVPMS']='azure'
#Morr_Processes_signed_colors['QMULTS']='cyan'
#Morr_Processes_signed_colors['QMULTR']='green'
#Morr_Processes_signed_colors['PRACS']='darkgreen'
#Morr_Processes_signed_colors['PSACWG']='palegreen'
#Morr_Processes_signed_colors['PGSACW']='gray'
#Morr_Processes_signed_colors['PGRACS']='orange'
#Morr_Processes_signed_colors['PRDG']='springgreen'
#Morr_Processes_signed_colors['EPRDG']='coral'
#Morr_Processes_signed_colors['EVPMG']='sage'
#Morr_Processes_signed_colors['PGMLT']='mediumpurple'
#Morr_Processes_signed_colors['PRACI']='lightsteelblue'
#Morr_Processes_signed_colors['PIACRS']='darkslategrey'
#Morr_Processes_signed_colors['PRACIS']='skyblue'
#Morr_Processes_signed_colors['EPRD']='beige'
#Morr_Processes_signed_colors['EPRDS']='saddlebrown'
#Morr_Processes_signed_colors['PRACG']='violet'
#Morr_Processes_signed_colors['QMULTG']='pink'
#Morr_Processes_signed_colors['QMULTRG']='indigo'
#Morr_Processes_signed_colors['MNUCCR']='lightcyan'
#Morr_Processes_signed_colors['MNUCCC']='indigo'
#Morr_Processes_signed_colors['MNUCCD']='indigo'
#Morr_Processes_signed_colors['PRAI']='lime'
#Morr_Processes_signed_colors['PRCI']='peru'
#Morr_Processes_signed_colors['PSACWS']='maroon'
#Morr_Processes_signed_colors['PIACR']='black'
#Morr_Processes_signed_colors['PSACWI']='gold'
#Morr_Processes_signed_colors['PSACR']='lightgray'
#Morr_Processes_signed_colors['EPCC']='lightblue'
#Morr_Processes_signed_colors['EPRCI']='blue'
#
#
##
#
#
#def load_wrf_morr_mass_proc(filename,add_coordinates=None,quantity='mixing ratio',slice_time=slice(None)):
#    from wrfload import loadwrfcube, derivewrfcube
#    from iris.cube import CubeList
#    Proclist=Proclist_Morr_mass
#    #Dict={}
#    cubelist_out=CubeList()
#    print('start loading processes')
#
#    if quantity=='volume':
#        print('start calculating density')
#        rho=derivewrfcube(filename,'density',slice_time=slice_time)
#    for i_process,process in enumerate(Proclist):
#        print(process)
#        if (i_process==0):
#            cube=loadwrfcube(filename,process+'3D',add_coordinates=add_coordinates,slice_time=slice_time)
#            #print(cube)
#            #print(add_coordinates)
#            if add_coordinates=='pz':
#                z_coord=cube.coord('geopotential')
#                p_coord=cube.coord('pressure')
#            if quantity=='volume':
#                cube=cube*rho
#            cube.rename(process)
#
#            #Dict[process]=cube
#            cubelist_out.append(cube)
#        else:
#             cube=loadwrfcube(filename,process+'3D',slice_time=slice_time)
#             if add_coordinates=='pz':
#                 cube.add_aux_coord(z_coord,(0,1,2,3))
#                 cube.add_aux_coord(p_coord,(0,1,2,3))
#             #Cubelist.append(cube)
#             if quantity=='volume':
#                 cube=cube*rho
#             cube.rename(process)
#             print(cube)
#             cubelist_out.append(cube)
#	       #Dict[process]=cube
#
#    return cubelist_out
#    #return Dict
#
#def load_wrf_morr_mass_processes_cubelist(filename,add_coordinates=None,quantity='mixing ratio',slice_time=slice(None)):
#    from wrfload import loadwrfcube, derivewrfcube
#    from iris import cube
#    Proclist=Proclist_Morr_mass
#    cubelist=cube.CubeList
#    print('start loading processes')
#
#    if quantity=='volume':
#        print('start calculating density')
#        rho=derivewrfcube(filename,'density',slice_time=slice_time)
#    for i_process,process in enumerate(Proclist):
#        print(process)
#        if (i_process==0):
#            cube=loadwrfcube(filename,process+'3D',add_coordinates=add_coordinates,slice_time=slice_time)
#            #print(cube)
#            #print(add_coordinates)
#            if add_coordinates=='pz':
#                z_coord=cube.coord('geopotential')
#                p_coord=cube.coord('pressure')
#            if quantity=='volume':
#                cube=cube*rho
#            cube.rename(process)
#
#            cubelist.append(cube)
#
#        else:
#             cube=loadwrfcube(filename,process+'3D',slice_time=slice_time)
#             if add_coordinates=='pz':
#                 cube.add_aux_coord(z_coord,(0,1,2,3))
#                 cube.add_aux_coord(p_coord,(0,1,2,3))
#             #Cubelist.append(cube)
#             if quantity=='volume':
#                 cube=cube*rho
#             cube.rename(process)
#             cubelist.append(cube)
#
#    return cubelist


#
#
#def load_wrf_morr_mass_proc_signed(filename,add_coordinates=None,quantity='mixing ratio',slice_time=slice(None),absolute_value=False):
#    from wrfload import loadwrfcube, derivewrfcube
#    from iris.cube import CubeList
#    Proclist=Proclist_Morr_mass
#    #Dict={}
#    cubelist_out=CubeList()
#    print('start loading processes')
#
#    List_signed=list(Morr_Processes_signed.keys())
#
#    if quantity=='volume':
#        print('start calculating density')
#        rho=derivewrfcube(filename,'density',slice_time=slice_time,add_coordinates=add_coordinates)
#        
#        
#        
#    for i_process,process in enumerate(Proclist):
#        print(process)
#        if (i_process==0):
#            cube=loadwrfcube(filename,process+'3D',add_coordinates=add_coordinates,slice_time=slice_time)
#            #print(cube)
#            #print(add_coordinates)
#            if add_coordinates=='pz':
#                z_coord=cube.coord('geopotential')
#                p_coord=cube.coord('pressure')
#            if quantity=='volume':
#                cube=cube*rho
#            cube.rename(process)
#
#            #Dict[process]=cube
#            cubelist_out.append(cube)
#
#        else:
#             if process in List_signed:
#                #Dict_1=split_sign_process(filename,process+'3D',Morr_Processes_signed[process][0],Morr_Processes_signed[process][1],add_coordinates=add_coordinates)
#                List_1=split_sign_process(filename,process+'3D',Morr_Processes_signed[process][0],Morr_Processes_signed[process][1],add_coordinates=add_coordinates)
#                if add_coordinates=='pz':
##                   for proc in Dict_1:
##                       Dict_1[proc].add_aux_coord(z_coord,(0,1,2,3))
##                       Dict_1[proc].add_aux_coord(p_coord,(0,1,2,3))
#                    for proc in List_1:
#                       proc.add_aux_coord(z_coord,(0,1,2,3))
#                       proc.add_aux_coord(p_coord,(0,1,2,3))
#
#                if quantity=='volume':
#                    #for proc in Dict_1:
#                    #    Dict_1[proc]= Dict_1[proc]*rho
#                    for proc in List_1:
#                        proc= proc*rho
#
#                cubelist_out.extend(List_1)
#             else:
#                cube=loadwrfcube(filename,process+'3D',slice_time=slice_time,add_coordinates=add_coordinates)
#                if add_coordinates=='pz':
#                    cube.add_aux_coord(z_coord,(0,1,2,3))
#                    cube.add_aux_coord(p_coord,(0,1,2,3))
#                if quantity=='volume':
#                    cube=cube*rho
#                cube.rename(process)
#                cubelist_out.append(cube)
#                #Dict[process]=cube
#
#    return cubelist_out
#    #return Dict

def load_wrf_variables_signed(filename,variable_list,split_dict,add_coordinates=None,constraint=None,quantity='mixing ratio',slice_time=slice(None),absolute_value=False,parallel_pool=None,debug_nproc=None):
    from wrfload import loadwrfcube, derivewrfcube
    from iris.cube import CubeList
    import numpy as np
    cubelist_out=CubeList()
    
    if (debug_nproc is not None):
        variable_list=variable_list[1:debug_nproc+1]

    List_signed=list(split_dict.keys())

    if quantity=='volume':
        rho=derivewrfcube(filename,'density',slice_time=slice_time,add_coordinates=add_coordinates,constraint=constraint)
    if add_coordinates=='pz':   
        cube=loadwrfcube(filename,variable_list[1],add_coordinates=add_coordinates,slice_time=slice_time,constraint=constraint)
        z_coord=cube.coord('geopotential')
        p_coord=cube.coord('pressure')
        
    for variable in variable_list:

        if variable in List_signed:
            List_1=split_sign_variable(filename,variable,split_dict[variable][0],split_dict[variable][1],add_coordinates=add_coordinates,constraint=constraint)
            if add_coordinates=='pz':
                for proc in List_1:
                   proc.add_aux_coord(z_coord,(0,1,2,3))
                   proc.add_aux_coord(p_coord,(0,1,2,3))
            if quantity=='volume':
                for i,variable in enumerate(List_1):
                    name=variable.name()
                    List_1[i]=(rho*variable)
                    List_1[i].rename(name) 
            cubelist_out.extend(List_1)
        else:
            cube=loadwrfcube(filename,variable,slice_time=slice_time,add_coordinates=add_coordinates,constraint=constraint)
            cube.data=np.abs(cube.data)
            if add_coordinates=='pz':
                cube.add_aux_coord(z_coord,(0,1,2,3))
                cube.add_aux_coord(p_coord,(0,1,2,3))
            if quantity=='volume':
                cube=cube*rho
            cube.rename(variable)
            cubelist_out.append(cube)
    return cubelist_out
#    
#    
#def load_wrf_proc_signed_morr(filename,process_list,split_dict,add_coordinates=None,constraint=None,quantity='mixing ratio',slice_time=slice(None),absolute_value=False,parallel_pool=None,debug_nproc=None):
#    from wrfload import loadwrfcube, derivewrfcube
#    from iris.cube import CubeList
#    #Dict={}
#    cubelist_out=CubeList()
#    #print('start loading processes')
#    
#
#    
#    if (debug_nproc is not None):
#        process_list=process_list[1:debug_nproc+1]
#
#    List_signed=list(split_dict.keys())
#
#    if quantity=='volume':
#        #print('start calculating density')
#        rho=derivewrfcube(filename,'density',slice_time=slice_time,add_coordinates=add_coordinates,constraint=constraint)
#    if add_coordinates=='pz':   
#        cube=loadwrfcube(filename,process_list[1],add_coordinates=add_coordinates,slice_time=slice_time,constraint=constraint)
#        z_coord=cube.coord('geopotential')
#        p_coord=cube.coord('pressure')
#        
#    for process in process_list:
#        if process in List_signed:
#            #Dict_1=split_sign_process(filename,process+'3D',Morr_Processes_signed[process][0],Morr_Processes_signed[process][1],add_coordinates=add_coordinates)
#            List_1=split_sign_process(filename,process,split_dict[process][0],split_dict[process][1],add_coordinates=add_coordinates,constraint=constraint)
#            if add_coordinates=='pz':
#                for proc in List_1:
#                   proc.add_aux_coord(z_coord,(0,1,2,3))
#                   proc.add_aux_coord(p_coord,(0,1,2,3))
#            if quantity=='volume':
#                for i,process in enumerate(List_1):
#                    name=process.name()
#                    List_1[i]=(rho*process)
#                    List_1[i].rename(name)
#                
#            cubelist_out.extend(List_1)
#        else:
#            cube=loadwrfcube(filename,process+'3D',slice_time=slice_time,add_coordinates=add_coordinates,constraint=constraint)
#            if add_coordinates=='pz':
#                cube.add_aux_coord(z_coord,(0,1,2,3))
#                cube.add_aux_coord(p_coord,(0,1,2,3))
#            if quantity=='volume':
#                cube=cube*rho
#            cube.rename(process)
#            #print(cube)
#            cubelist_out.append(cube)
#
#    return cubelist_out
#
#def load_wrf_morr_mass_proc_individual(filename,process,add_coordinates=None,quantity='mixing ratio',slice_time=slice(None)):
#    from wrfload import loadwrfcube, derivewrfcube
#    Dict={}
#    #if quantity=='volume':
#    #    rho=derivewrfcube(filename,'density',slice_time=slice_time)
#
##    if process=='PCC':
##        cube=loadwrfcube(filename,process+'3D',add_coordinates=add_coordinates,slice_time=slice_time)
##        cube1=cube.copy()
##        cube2=cube.copy()
##        cube1.data=np.clip(cube1.data,a_min=-np.inf,a_max=0)
##        if quantity=='volume':
##            cube1=cube1*rho
##        cube1.rename('PCC')
##        Dict['PCC']=cube1
##        cube2=cube[:]
##        cube2.data=np.clip(cube.data,a_min=0,a_max=np.inf)
##        if quantity=='volume':
##            cube2=cube2*rho
##        cube2.rename('EPCC')
##        Dict['EPCC']=cube2
##
##    else:
#    cube=loadwrfcube(filename,process+'3D',add_coordinates=add_coordinates,slice_time=slice_time)
##    if add_coordinates=='pz':
##        cube.add_aux_coord(z_coord,(0,1,2,3))
##        cube.add_aux_coord(p_coord,(0,1,2,3))
##  #Cubelist.append(cube)
##    if quantity=='volume':
##        cube=cube*rho
#    cube.rename(process)
#
#    Dict[process]=cube
#
#    #return Cubelist
#    return Dict




#def load_wrf_morr_num_proc(filename,add_coordinates=None,quantity='volume'):
#    from wrfload import loadwrfcube, derivewrfcube
#    Dict={}
#    if quantity=='volume':
#        rho=derivewrfcube(filename,'density')
#
#
#    #Cubelist=[]
#    Dict={}
#    for i_process,process in enumerate(Proclist_Morr_number):
#        if (i_process==0):
#            cube=loadwrfcube(filename,process+'3D',add_coordinates=add_coordinates)
#            cube.rename(process)
#            #Cubelist.append(cube)
#            if add_coordinates=='pz':
#                z_coord=cube.coord('geopotential')
#                p_coord=cube.coord('pressure')
#        else:
#            if process=='NSMLTR':
#                cube=loadwrfcube(filename,process)
#            else:
#                cube=loadwrfcube(filename,process+'3D')
#            cube.rename(process)
#            if add_coordinates=='pz':
#                cube.add_aux_coord(z_coord,(0,1,2,3))
#                cube.add_aux_coord(p_coord,(0,1,2,3))
#
#        #Cubelist.append(cube)
#        if quantity=='volume':
#            cube=cube*rho
#
#        Dict[process]=cube
    #return Cubelist
#    return Dict
#
#Hydropath_list=[
#    'VAPORCLOUD',
#    'VAPORRAIN',
#    'VAPORICE',
#    'VAPORSNOW',
#    'VAPORGRAUP',
#    'CLOUDRAIN',
#    'CLOUDICE',
#    'CLOUDSNOW',
#    'CLOUDGRAUP',
#    'RAINICE',
#    'RAINSNOW',
#    'RAINGRAUP',
#    'ICESNOW',
#    'ICEGRAUP',
#    'SNOWGRAUP']
#
#def calculate_wrf_morr_path_hydrometeors(filename,slice_time=slice(None)):
#    Dict={}
#    #Cubelist=[]
#    for path in Hydropath_list:
#        cube=calculate_wrf_morr_path(filename,path,slice_time=slice_time)
#        #Cubelist.append(cube)
#        Dict[path]=cube
#    #return Cubelist
#    return Dict
#
#Phasepath_list=['vaporliquid','vaporfrozen','liquidfrozen']
#
#def calculate_wrf_morr_path_phases(filename,slice_time=slice(None)):
#    Dict={}
#    #Cubelist=[]
#    for path in Phasepath_list:
#        print('loading ',  path)
#        cube=calculate_wrf_morr_path(filename,path,slice_time=slice_time)
#        print(path, ' loaded')
#
#        #Cubelist.append(cube)
#        Dict[path]=cube
#        print(Dict[path].data)
#    #return Cubelist
#    return Dict

def calculate_wrf_thompson_path(filename,path,add_coordinates=None,quantity='volume'):
    if (path=='processes_mass'):
        out=load_wrf_thom_mass_proc(filename,add_coordinates)
    if (path=='processes_number'):
        out=load_wrf_thom_number_proc(filename,add_coordinates)
    else:
        print('option not avaliable')
    return out


thompson_processes_mass= list(List_Processes_Thompson_Mass)
#thompson_processes_mass_split={}
thompson_processes_number= list(List_Processes_Thompson_Number)
thompson_processes_number_split={}

def calculate_wrf_mp_path(filename,processes=None,microphysics_scheme=None, signed=False,constraint=None,add_coordinates=None,quantity='volume',slice_time=slice(None),parallel_pool=None,debug_nproc=None):

    if microphysics_scheme=='morrison':
        if processes=='mass':
            process_list=morrison_processes_mass
            #process_list.remove('PCCN')
            if signed==True:
                split_dict=morrison_processes_mass_split
            else:
                split_dict={}
        elif processes=='number':
            process_list=morrison_processes_number
            if signed==True:
                split_dict='morrison_processes_number_split'
            else:
                split_dict={}    
#                
#    variablelist=variable_list(filename)
#    if 'PCC3D' in variablelist:
#        add_str='3D'
#        
#        for i,process in enumerate(process_list)
#            process_list[i]=process+add_str
#        
#        for i,process in enumerate(split_dict)
#            split_dict[i]=process+add_str


        cube_list_out=load_wrf_variables_signed(filename,variable_list=process_list,split_dict=split_dict,add_coordinates=add_coordinates,quantity=quantity,slice_time=slice_time,constraint=constraint,parallel_pool=parallel_pool,debug_nproc=debug_nproc)

    elif microphysics_scheme=='thompson':
        if processes=='mass':
            process_list=thompson_processes_mass
            if signed==True:
                split_dict=thompson_processes_mass_split
            else:
                split_dict={}
        elif processes=='number':
            process_list=thompson_processes_number
            if signed==True:
                split_dict=thompson_processes_number_split
            else:
                split_dict={}
        cube_list_out=load_wrf_variables_signed(filename,variable_list=process_list,split_dict=split_dict,add_coordinates=add_coordinates,quantity=quantity,slice_time=slice_time,constraint=constraint,parallel_pool=parallel_pool,debug_nproc=debug_nproc)
 
    return cube_list_out
    

    
#def calculate_wrf_morr_path(filename,path,add_coordinates=None,quantity='volume',slice_time=slice(None)):
#    if (path=='processes_mass'):
#        out=load_wrf_morr_mass_proc(filename,add_coordinates,quantity,slice_time=slice_time)
#    if (path=='processes_number'):
#        out=load_wrf_morr_num_proc(filename,add_coordinates,quantity)
#    if path=='hydrometeor':
#        out=calculate_wrf_morr_path_hydrometeors(filename,slice_time=slice_time)
#    if path=='phase':
#        out=calculate_wrf_morr_path_phases(filename,slice_time=slice_time)
#    if (path=='vaporliquid'):
#        out=calculate_wrf_morr_path_vaporliquid(filename,slice_time=slice_time)
#    if path=='vaporfrozen':
#        out=calculate_wrf_morr_path_vaporfrozen(filename,slice_time=slice_time)
#    if path=='liquidfrozen':
#        out=calculate_wrf_morr_path_liquidfrozen(filename,slice_time=slice_time)
#    if path=='VAPORCLOUD':
#        out=calculate_wrf_morr_path_VAPORCLOUD(filename,slice_time=slice_time)
#    if path=='VAPORRAIN':
#        out=calculate_wrf_morr_path_VAPORRAIN(filename,slice_time=slice_time)
#    if path=='VAPORICE':
#        out=calculate_wrf_morr_path_VAPORICE(filename,slice_time=slice_time)
#    if path=='VAPORSNOW':
#        out=calculate_wrf_morr_path_VAPORSNOW(filename,slice_time=slice_time)
#    if path=='VAPORGRAUP':
#        out=calculate_wrf_morr_path_VAPORGRAUP(filename,slice_time=slice_time)
#    if path=='CLOUDRAIN':
#        out=calculate_wrf_morr_path_CLOUDRAIN(filename,slice_time=slice_time)
#    if path=='CLOUDICE':
#        out=calculate_wrf_morr_path_CLOUDICE(filename,slice_time=slice_time)
#    if path=='CLOUDSNOW':
#        out=calculate_wrf_morr_path_CLOUDSNOW(filename,slice_time=slice_time)
#    if path=='CLOUDGRAUP':
#        out=calculate_wrf_morr_path_CLOUDGRAUP(filename,slice_time=slice_time)
#    if path=='RAINICE':
#        out=calculate_wrf_morr_path_RAINICE(filename,slice_time=slice_time)
#    if path=='RAINSNOW':
#        out=calculate_wrf_morr_path_RAINSNOW(filename,slice_time=slice_time)
#    if path=='RAINGRAUP':
#        out=calculate_wrf_morr_path_RAINGRAUP(filename,slice_time=slice_time)
#    if path=='ICESNOW':
#        out=calculate_wrf_morr_path_ICESNOW(filename,slice_time=slice_time)
#    if path=='ICEGRAUP':
#        out=calculate_wrf_morr_path_ICEGRAUP(filename,slice_time=slice_time)
#    if path=='SNOWGRAUP':
#        out=calculate_wrf_morr_path_SNOWGRAUP(filename,slice_time=slice_time)
#    else:
#        print('path string unknown')
#    return out
#
#def calculate_wrf_morr_latentheating(filename):
#    #Load and add up all process rates between water vapour and cloud droplets
#    from wrfload import loadwrfcube
#    LHREVP=loadwrfcube(filename,'LHREVP')
#    LHRFRZ=loadwrfcube(filename,'LHRFRZ')
#    LHRSUB=loadwrfcube(filename,'LHRSUB')
#    LHR=-1*LHREVP+LHRFRZ+LHRSUB
#    LHR.rename('latent heating rate')
#    return LHR
#
#
#
#def calculate_wrf_morr_path_VAPORCLOUD(filename,slice_time=slice(None)):
#    #Load and add up all process rates between water vapour and cloud droplets
#    from wrfload import loadwrfcube
#    #print('calculate process rates VAPOR/CLOUD')
#    PCC= loadwrfcube(filename, 'PCC3D',slice_time=slice_time)
#    #PCCN=loadwrfcube(filename, 'PCCN3D')
#    P_VAPORCLOUD = PCC#+PCCN
#    P_VAPORCLOUD.rename('P_VAPORCLOUD')
#    return P_VAPORCLOUD
#
#def calculate_wrf_morr_path_VAPORRAIN(filename,slice_time=slice(None)):
#    #Load and add up all process rates between water vapour and cloud droplets
#    from wrfload import loadwrfcube
#    #print('calculate process rates VAPOR/RAIN')
#    P_VAPORRAIN = loadwrfcube(filename, 'PRE3D',slice_time=slice_time)  #EVAP OF RAIN
#    P_VAPORRAIN.rename('P_VAPORRAIN')
#
#    return P_VAPORRAIN
#
#def calculate_wrf_morr_path_VAPORICE(filename,slice_time=slice(None)):
#    #Load and add up all process rates between water vapour and cloud droplets
#    from wrfload import loadwrfcube
#    #print('calculate process rates VAPOR/ICE')
#    PRD=loadwrfcube(filename, 'PRD3D',slice_time=slice_time)     # DEP CLOUD ICE
#    EPRD=loadwrfcube(filename, 'EPRD3D',slice_time=slice_time)     # SUBLIMATION CLOUD ICE
#    MNUCCD=loadwrfcube(filename, 'MNUCCD3D',slice_time=slice_time) # CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
#    P_VAPORICE = PRD + EPRD + MNUCCD
#    P_VAPORICE.rename('P_VAPORICE')
#    return P_VAPORICE
#
#def calculate_wrf_morr_path_VAPORSNOW(filename,slice_time=slice(None)):
#    #Load and add up all process rates between water vapour and cloud droplets
#    from wrfload import loadwrfcube
#    #print('calculate process rates VAPOR/SNOW')
#    EVPMS=loadwrfcube(filename, 'EVPMS3D',slice_time=slice_time) # CHANGE Q MELTING SNOW EVAPORATING
#    EPRDS=loadwrfcube(filename, 'EPRDS3D',slice_time=slice_time) #    SUBLIMATION SNOW
#    P_VAPORSNOW=EVPMS+EPRDS
#    P_VAPORSNOW.rename('P_VAPORSNOW')
#    return P_VAPORSNOW
#
#def calculate_wrf_morr_path_VAPORGRAUP(filename,slice_time=slice(None)):
#    #Load and add up all process rates between water vapour and cloud droplets
#    from wrfload import loadwrfcube
#    #print('calculate process rates VAPOR/GRAUPEL')
#    EVPMG= loadwrfcube(filename, 'EVPMG3D',slice_time=slice_time)  # CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
#    PRDG=loadwrfcube(filename,'PRDG3D',slice_time=slice_time)  # DEP OF GRAUPEL
#    EPRDG=loadwrfcube(filename, 'EPRDG3D',slice_time=slice_time)  #  SUB OF GRAUPEL
#    P_VAPORGRAUPEL=EVPMG+PRDG+EPRDG
#    P_VAPORGRAUPEL.rename('P_VAPORGRAUPEL')
#    return P_VAPORGRAUPEL
#
#def calculate_wrf_morr_path_CLOUDRAIN(filename,slice_time=slice(None)):
#    #Load and add up all process rates between cloud droplets and rain
#    from wrfload import loadwrfcube
#    #print('calculate process rates CLOUD/RAIN')
#    PRA= loadwrfcube(filename, 'PRA3D',slice_time=slice_time)      # ACCRETION DROPLETS BY RAIN
#    PRC= loadwrfcube(filename, 'PRC3D',slice_time=slice_time)    # AUTOCONVERSION DROPLETS
#    P_CLOUDRAIN=PRA+PRC
#    P_CLOUDRAIN.rename('P_CLOUDRAIN')
#    return P_CLOUDRAIN
#
#
#def calculate_wrf_morr_path_CLOUDICE(filename,slice_time=slice(None)):
#    #Load and add up all process rates between ckoud droplets and cloud ice
#    from wrfload import loadwrfcube
#    # print('calculate process rates CLOUD/ICE')
#    PSACWI=loadwrfcube(filename, 'PSACWI3D',slice_time=slice_time)   # CHANGE Q DROPLET ACCRETION BY CLOUD ICE
#    QMULTS=loadwrfcube(filename, 'QMULTS3D',slice_time=slice_time)  # CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
#    QMULTG=loadwrfcube(filename, 'QMULTG3D',slice_time=slice_time)   # CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
#    P_CLOUDICE =PSACWI+QMULTS+QMULTG
#    P_CLOUDICE.rename('P_CLOUDICE')
#    return P_CLOUDICE
#
#
#def calculate_wrf_morr_path_CLOUDSNOW(filename,slice_time=slice(None)):
#    #Load and add up all process rates between cloud droplets and snow
#    from wrfload import loadwrfcube
#    # print('calculate process rates CLOUD/ICE')
#    PSACWS= loadwrfcube(filename, 'PSACWS3D',slice_time=slice_time)      # CHANGE Q DROPLET ACCRETION BY SNOW
#    P_CLOUDSNOW=PSACWS
#    P_CLOUDSNOW.rename('P_CLOUDSNOW')
#    return P_CLOUDSNOW
#
#def calculate_wrf_morr_path_CLOUDGRAUP(filename,slice_time=slice(None)):
#    #Load and add up all process rates between cloud droplets and graupel
#    from wrfload import loadwrfcube
#    # print('calculate process rates CLOUD/GRAUP')
#    PSACWG=loadwrfcube(filename, 'PSACWG3D',slice_time=slice_time)      #  CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
#    PGSACW=loadwrfcube(filename, 'PGSACW3D',slice_time=slice_time)   # CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
#    P_CLOUDGRAUP = PGSACW + PSACWG
#    P_CLOUDGRAUP.rename('P_CLOUDGRAUP')
#    return P_CLOUDGRAUP
#
#
#def calculate_wrf_morr_path_RAINICE(filename,slice_time=slice(None)):
#    #Load and add up all process rates between rain and cloud ice
#    from wrfload import loadwrfcube
#    # print('calculate process rates RAIN/ICE')
#    QMULTR=loadwrfcube(filename, 'QMULTR3D',slice_time=slice_time)      # CHANGE Q DUE TO ICE RAIN/SNOW
#    QMULTRG=loadwrfcube(filename, 'QMULTRG3D',slice_time=slice_time)                         # CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL
#    P_RAINICE = QMULTR+QMULTRG
#    P_RAINICE.rename('P_RAINICE')
#    return P_RAINICE
#
#def calculate_wrf_morr_path_RAINSNOW(filename,slice_time=slice(None)):
#    #Load and add up all process rates between rain and snow
#    from wrfload import loadwrfcube
#    # print('calculate process rates RAIN/SNOW')
#    PSMLT=loadwrfcube(filename, 'PSMLT3D',slice_time=slice_time)     # CHANGE Q MELTING SNOW TO RAIN
#    PIACRS=loadwrfcube(filename, 'PIACRS3D',slice_time=slice_time)                           #CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
#    P_RAINSNOW =PSMLT+PIACRS
#    P_RAINSNOW.rename('P_RAINSNOW')
#    return P_RAINSNOW
#
#def calculate_wrf_morr_path_RAINGRAUP(filename,slice_time=slice(None)):
#    #Load and add up all process rates between rain and graupel
#    from wrfload import loadwrfcube
#    # print('calculate process rates RAIN/GRAUPEL')
#    MNUCCR=loadwrfcube(filename, 'MNUCCR3D',slice_time=slice_time)     # CHANGE Q DUE TO CONTACT FREEZ RAIN
#    PIACR=loadwrfcube(filename, 'PIACR3D',slice_time=slice_time)      # CHANGE QR, ICE-RAIN COLLECTION
#    PRACG=loadwrfcube(filename, 'PRACG3D',slice_time=slice_time)      #CHANGE IN Q COLLECTION RAIN BY GRAUPEL
#    PGRACS=loadwrfcube(filename, 'PGRACS3D',slice_time=slice_time) # CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
#    PGMLT=loadwrfcube(filename, 'PGMLT3D',slice_time=slice_time) #  CHANGE Q MELTING OF GRAUPEL
#    P_RAINGRAUP =MNUCCR+PIACR+PRACG+PGRACS+PGMLT
#    P_RAINGRAUP.rename('P_RAINGRAUP')
#    return P_RAINGRAUP
#
#def calculate_wrf_morr_path_ICESNOW(filename,slice_time=slice(None)):
#    #Load and add up all process rates between cloud ice and snow
#    from wrfload import loadwrfcube
#    # print('calculate process rates ICE/SNOW')
#    PRAI = loadwrfcube(filename, 'PRAI3D',slice_time=slice_time)      # CHANGE Q ACCRETION CLOUD ICE BY SNOW
#    PRCI=loadwrfcube(filename, 'PRCI3D',slice_time=slice_time)      # CHANGE Q AUTOCONVERSIN CLOUD ICE TO SNOW
#    PRACIS=loadwrfcube(filename, 'PRACIS3D',slice_time=slice_time)     # CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
#    P_ICESNOW = PRAI + PRCI + PRACIS
#    P_ICESNOW.rename('P_ICESNOW')
#    return P_ICESNOW
#
#def calculate_wrf_morr_path_ICEGRAUP(filename,slice_time=slice(None)):
#    #Load and add up all process rates between cloud ice and graupel
#    from wrfload import loadwrfcube
#    # print('calculate process rates ICE/GRAUPEL')
#    PRACI=loadwrfcube(filename, 'PRACI3D',slice_time=slice_time)     # CHANGE QI, ICE-RAIN COLLECTION
#    P_ICEGRAUP = PRACI
#    P_ICEGRAUP .rename('P_ICEGRAUP ')
#    return P_ICEGRAUP
#
#def calculate_wrf_morr_path_SNOWGRAUP(filename,slice_time=slice(None)):
#    #Load and add up all process rates between snow and graupel
#    from wrfload import loadwrfcube
#    # print('calculate process rates SNOW/GRAUPEL')
#    P_SNOWGRAUP = 0*loadwrfcube(filename, 'PRACI3D',slice_time=slice_time)   # Dummy zeros, since no pathway process found yet
#    P_SNOWGRAUP.rename('P_SNOWGRAUP')
#    return P_SNOWGRAUP
#
#
#def calculate_wrf_morr_path_vaporliquid(filename,slice_time=slice(None)):
#    #Load and add up all process rates between ice phase and water vapour:
#    #print('calculate processes deposition/sublimation')
#    PVAPORCLOUD=calculate_wrf_morr_path_VAPORCLOUD(filename,slice_time=slice_time)
#    PVAPORRAIN=calculate_wrf_morr_path_VAPORRAIN(filename,slice_time=slice_time)
#    P_vaporliquid=  PVAPORCLOUD + PVAPORRAIN
#    P_vaporliquid.rename('PVAPORLIQUID')
#    return P_vaporliquid
#
#def calculate_wrf_morr_path_vaporfrozen(filename,slice_time=slice(None)):
#    #Load and add up all process rates between ice phase and water vapour:
#    #print('calculate processes deposition/sublimation')
#    PVAPORICE=calculate_wrf_morr_path_VAPORICE(filename,slice_time=slice_time)
#    PVAPORSNOW=calculate_wrf_morr_path_VAPORSNOW(filename,slice_time=slice_time)
#    PVAPORGRAUP=calculate_wrf_morr_path_VAPORGRAUP(filename,slice_time=slice_time)
#    P_vaporfrozen=PVAPORICE+PVAPORSNOW+PVAPORGRAUP
#    P_vaporfrozen.rename('PVAPORFROZEN')
#    return P_vaporfrozen
#
#def calculate_wrf_morr_path_liquidfrozen(filename,slice_time=slice(None)):
#    #Load and add up all process rates between frozen and liquid phase
#    #print('calculate processes freezing/melting')
#    PCLOUDICE=calculate_wrf_morr_path_CLOUDICE(filename,slice_time=slice_time)
#    PRAINICE=calculate_wrf_morr_path_RAINICE(filename,slice_time=slice_time)
#    PCLOUDSNOW=calculate_wrf_morr_path_CLOUDSNOW(filename,slice_time=slice_time)
#    PRAINSNOW=calculate_wrf_morr_path_RAINSNOW(filename,slice_time=slice_time)
#    PRAINGRAUP=calculate_wrf_morr_path_RAINGRAUP(filename,slice_time=slice_time)
#    PCLOUDGRAUP=calculate_wrf_morr_path_CLOUDGRAUP(filename,slice_time=slice_time)
#    P_liquidfrozen=PCLOUDICE+PRAINICE+PCLOUDSNOW+PRAINSNOW+PRAINGRAUP+PCLOUDGRAUP
#    P_liquidfrozen.rename('PLIQUIDFROZEN')
#
#    return P_liquidfrozen
#
#def sum_cubes(filename,name,list_names):
#    from wrfload import loadwrfcube
#    P_out=loadwrfcube(filename, list_names[0])
#    for name_i in list_names[1:]:
#        P_out=P_out+loadwrfcube(filename, name_i)
#    P_out.rename(name)
#    return P_out
#
#
#
#
#def calculate_wrf_thom_path_VAPORCLOUD(filename):
#    list_names=['']
#    name='P_VAPORCLOUD'
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_VAPORRAIN(filename):
#    #Load and add up all process rates between water vapour and cloud droplets
#    list_names=['']
#    name='P_VAPORRAIN'
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_VAPORICE(filename):
#    #Load and add up all process rates between water vapour and cloud droplets
#    list_names=['']
#    name='P_VAPORICE'
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_VAPORSNOW(filename):
#    list_names=['']
#    name='P_VAPORSNOW'
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_VAPORGRAUP(filename):
#    #Load and add up all process rates between water vapour and cloud droplets
#    list_names=['']
#    name='P_VAPORGRAUPEL'
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_CLOUDRAIN(filename):
#     #Load and add up all process rates between water vapour and cloud droplets
#    list_names=['']
#    name='P_CLOUDRAIN'
#    return sum_cubes(filename,name,list_names)
#
#
#def calculate_wrf_thom_path_CLOUDICE(filename):
#    #Load and add up all process rates between ckoud droplets and cloud ice
#    list_names=['']
#    name='P_CLOUDICE'
#    return sum_cubes(filename,name,list_names)
#
#
#
#def calculate_wrf_thom_path_CLOUDSNOW(filename):
#    #Load and add up all process rates between cloud droplets and snow
#    list_names=['']
#    name='P_CLOUDSNOW'
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_CLOUDGRAUP(filename):
#    #Load and add up all process rates between cloud droplets and graupel
#    list_names=['']
#    name='P_CLOUDGRAUP'
#    return sum_cubes(filename,name,list_names)
#
#
#def calculate_wrf_thom_path_RAINICE(filename):
#    #Load and add up all process rates between rain and cloud ice
#    list_names=['']
#    name='P_RAINICE'
#    return sum_cubes(filename,name,list_names)
#def calculate_wrf_thom_path_RAINSNOW(filename):
#    #Load and add up all process rates between rain and snow
#    list_names=['']
#    name='P_RAINICE'
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_RAINGRAUP(filename):
#    #Load and add up all process rates between rain and graupel
#    list_names=['']
#    name='P_RAINGRAUP'
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_ICESNOW(filename):
#    #Load and add up all process rates between cloud ice and snow
#    list_names=['']
#    name='P_ICESNOW'
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_ICEGRAUP(filename):
#    #Load and add up all process rates between cloud ice and graupel
#    list_names=['']
#    name='P_ICEGRAUP'
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_SNOWGRAUP(filename):
#    #Load and add up all process rates between snow and graupel
#    list_names=['']
#    name='P_SNOWGRAUP'
#    return sum_cubes(filename,name,list_names)
#
#
#def calculate_wrf_thom_path_vaporliquid(filename):
#    #Load and add up all process rates between ice phase and water vapour:
#    #print('calculate processes deposition/sublimation')
#    PVAPORCLOUD=calculate_wrf_thom_path_VAPORCLOUD(filename)
#    PVAPORRAIN=calculate_wrf_thom_path_VAPORRAIN(filename)
#    name=('PVAPORLIQUID')
#    return sum_cubes(filename,name,list_names)
#
#def calculate_wrf_thom_path_vaporfrozen(filename):
#    #Load and add up all process rates between ice phase and water vapour:
#    #print('calculate processes deposition/sublimation')
#    PVAPORICE=calculate_wrf_thom_path_VAPORICE(filename)
#    PVAPORSNOW=calculate_wrf_thom_path_VAPORSNOW(filename)
#    PVAPORGRAUP=calculate_wrf_thom_path_VAPORGRAUP(filename)
#    P_vaporfrozen=PVAPORICE+PVAPORSNOW+PVAPORGRAUP
#    P_vaporfrozen.rename('PVAPORFROZEN')
#    return P_vaporfrozen
#
#def calculate_wrf_thom_path_liquidfrozen(filename):
#    #Load and add up all process rates between frozen and liquid phase
#    #print('calculate processes freezing/melting')
#    PCLOUDICE=calculate_wrf_thom_path_CLOUDICE(filename)
#    PRAINICE=calculate_wrf_thom_path_RAINICE(filename)
#    PCLOUDSNOW=calculate_wrf_thom_path_CLOUDSNOW(filename)
#    PRAINSNOW=calculate_wrf_thom_path_RAINSNOW(filename)
#    PRAINGRAUP=calculate_wrf_thom_path_RAINGRAUP(filename)
#    PCLOUDGRAUP=calculate_wrf_thom_path_CLOUDGRAUP(filename)
#    P_liquidfrozen=PCLOUDICE+PRAINICE+PCLOUDSNOW+PRAINSNOW+PRAINGRAUP+PCLOUDGRAUP
#    P_liquidfrozen.rename('PLIQUIDFROZEN')
#
#    return P_liquidfrozen
#
#

def processes_colors(microphysics_scheme=None,colors_processes='all'):
    Processes_signed_colors={}
    Processes_signed_names={}

    #for process in Processes:
        #Processes_signed_colors[process.name()]='gray'
  
        #set colors for specific processes: 
    if microphysics_scheme=='morrison':
        
        if colors_processes=='lumped':
            Processes_signed_colors=lumped_colors_morrison
        else:
         
            Processes_signed_colors['PCC']='lightblue'        
            Processes_signed_names['PCC']='Cond. droplets (PCC)'
    
            Processes_signed_colors['EPCC']='lightsalmon'
            Processes_signed_names['EPCC']='Evap droplets (EPCC)'
    
            Processes_signed_colors['PRE']='salmon'
            Processes_signed_names['PRE']='Evap. rain (PRE)'#''
            
            Processes_signed_colors['PRA']='lightgreen'
            Processes_signed_names['PRA']='Accretion (PRA)'#''
            
            Processes_signed_colors['PRC']='brown'         
            Processes_signed_names['PRC']='Autoconversion (PRC)'         
    
            Processes_signed_colors['PRDG']='cyan'
            Processes_signed_names['PRDG']='Dep. graupel (PRDG)'#''
    
            Processes_signed_colors['EPRDG']='gold'
            Processes_signed_names['EPRDG']='Subl. graupel (EPRDG)'#''
    
            Processes_signed_colors['EPRD']='khaki'
            Processes_signed_names['EPRD']='Subl. ice (EPRD)'# ''
    
            Processes_signed_colors['EPRDS']='wheat'
            Processes_signed_names['EPRDS']='Subl. snow (EPRDS)'#''
    
            Processes_signed_colors['PGRACS']='orange'
            Processes_signed_names['PGRACS']='Coll. rain/snow (PGRACS)'# ''
    
            Processes_signed_colors['PRDS']='darkorange'
            Processes_signed_names['PRDS']='Dep. snow (PRDS)'#'Dep. snow'
    
            Processes_signed_colors['PSACWG']='palegreen'
            Processes_signed_names['PSACWG']='Coll. droplets/graupel (PSACWG)'#'Coll. droplets/graupel'
    
    
    
            if colors_processes=='all':
                Processes_signed_colors['PRC']='brown'       
                Processes_signed_names['PRC']='PRC'         
    
                Processes_signed_colors['PRD']='lightseagreen'
                Processes_signed_names['PRD']='PRD'
    
                Processes_signed_colors['PCCN']='blue'
                Processes_signed_names['PCCN']='PCCN'
    
                Processes_signed_colors['PSMLT']='darkblue'
                Processes_signed_names['PSMLT']='PSMLT'
    
                Processes_signed_colors['EVPMS']='azure'
                Processes_signed_names['EVPMS']='EVPMS'
    
                Processes_signed_colors['QMULTS']='cyan'
                Processes_signed_names['QMULTS']='QMULTS'
    
                Processes_signed_colors['QMULTR']='green'
                Processes_signed_names['QMULTR']='QMULTR'
    
                Processes_signed_colors['PRACS']='darkgreen'
                Processes_signed_names['PRACS']='PRACS'
    
                Processes_signed_colors['PGSACW']='gray'
                Processes_signed_names['PGSACW']='PGSACW'
    
                Processes_signed_colors['PRDG']='springgreen'
                Processes_signed_names['PRDG']='PRDG'
    
                Processes_signed_colors['EPRDG']='coral'
                Processes_signed_names['EPRDG']='EPRDG'
    
                Processes_signed_colors['EVPMG']='mediumseagreen'
                Processes_signed_names['EVPMG']='EVPMG'
    
                Processes_signed_colors['PGMLT']='mediumpurple'
                Processes_signed_names['PGMLT']='PGMLT'
    
                Processes_signed_colors['PRACI']='lightsteelblue'
                Processes_signed_names['PRACI']='PRACI'
    
                Processes_signed_colors['PIACRS']='darkslategrey'
                Processes_signed_names['PIACRS']='PIACRS'
    
                Processes_signed_colors['PRACIS']='skyblue'
                Processes_signed_names['PRACIS']='PRACIS'
    
                Processes_signed_colors['EPRD']='beige'
                Processes_signed_names['EPRD']='EPRD'
    
                Processes_signed_colors['PRACG']='violet'
                Processes_signed_names['PRACG']='PRACG'
    
                Processes_signed_colors['QMULTG']='pink'
                Processes_signed_names['QMULTG']='QMULTG'
    
                Processes_signed_colors['QMULTRG']='indigo'
                Processes_signed_names['QMULTRG']='QMULTRG'
    
                Processes_signed_colors['MNUCCR']='lightcyan'
                Processes_signed_names['MNUCCR']='MNUCCR'
    
                Processes_signed_colors['MNUCCC']='indigo'
                Processes_signed_names['MNUCCC']='MNUCCC'
    
                Processes_signed_colors['MNUCCD']='indigo'
                Processes_signed_names['MNUCCD']='MNUCCD'
    
                Processes_signed_colors['PRAI']='lime'
                Processes_signed_names['PRAI']='PRAI'
    
                Processes_signed_colors['PRCI']='peru'
                Processes_signed_names['PRCI']='PRCI'
    
                Processes_signed_colors['PSACWS']='maroon'
                Processes_signed_names['PSACWS']='PSACWS'
    
                Processes_signed_colors['PIACR']='black'
                Processes_signed_names['PIACR']='PIACR'
    
                Processes_signed_colors['PSACWI']='gold'
                Processes_signed_names['PSACWI']='PSACWI'
    
                Processes_signed_colors['PSACR']='lightgray'
                Processes_signed_names['PSACR']='PSACR'
    
                Processes_signed_colors['EPRCI']='blue'
                Processes_signed_names['EPRCI']='EPRCI'
                
                Processes_signed_colors['QICF']='#99d8c9'
                Processes_signed_names['QICF']='QICF'

                Processes_signed_colors['QGRF']='#756bb1'
                Processes_signed_names['QGRF']='QGRF'

                Processes_signed_colors['QNIRF']='#bcbddc'
                Processes_signed_names['QNIRF']='QNIRF'
                
                Processes_signed_colors['QIIM']='#d95f0e'
                Processes_signed_names['QIIM']='QIIM'


    if microphysics_scheme=='thompson':
        
        
        if colors_processes=='lumped':
            Processes_signed_colors=lumped_colors_thompson
        else:
            
            Processes_signed_colors['PRI_INU']='lightseagreen'   
            Processes_signed_names['PRI_INU']='Ice nucleation (PRI_INU)'   
    
            Processes_signed_colors['PRR_RCW']='lightgreen'  
            Processes_signed_names['PRR_RCW']='Coll. droplets/rain (PRR_RCW)'  
    
            Processes_signed_colors['PRR_GML']='wheat'  
            Processes_signed_names['PRR_GML']='Melt. graupel (PRR_GML)' 
    
            Processes_signed_colors['PRI_IDE']='blue'   
            Processes_signed_names['PRI_IDE']='Dep. ice (PRI_IDE)'   
    
            Processes_signed_colors['PRS_SDE']='lightsteelblue' 
            Processes_signed_names['PRS_SDE']='Dep. snow (PRS_SDE)'
            
            Processes_signed_colors['PRS_IAU']='green'   
            Processes_signed_names['PRS_IAU']='Autoconv.ice (PRS_IAU)'   
    
            Processes_signed_colors['PRI_WFZ']='darkred'  
            Processes_signed_names['PRI_WFZ']='Freez. droplets (PRI_WFZ)'   
    
            Processes_signed_colors['PRG_GDE']='indigo'   
            Processes_signed_names['PRG_GDE']='Dep. graupel (PRG_GDE)'   
    
            Processes_signed_colors['E_PRS_SDE']='darkslateblue'  
            Processes_signed_names['E_PRS_SDE']='Subl. snow (E_PRS_SDE)'  
    
            Processes_signed_colors['PRG_RCG']='black'  
            Processes_signed_names['PRG_RCG']='Coll. rain/graupel (PRG_RCG)'
    
            Processes_signed_colors['PRV_REV']='lightsalmon'   #  Vapor->Water
            Processes_signed_names['PRV_REV']='Evaporation of rain (PRV_REV)'
    
            Processes_signed_colors['PRW_VCD']='lightblue'   #  Vapor->Water
            Processes_signed_names['PRW_VCD']='Condensation (PRW_VCD)'
    
            Processes_signed_colors['E_PRW_VCD']='coral'   #  Vapor->Water
            Processes_signed_names['E_PRW_VCD']='Evaporation (E_PRW_VCD)'
    
    
           
            if colors_processes=='all':        
                Processes_signed_colors['PRR_RCI']='magenta'   #  Vapor->Water
                Processes_signed_names['PRR_RCI']='PRR_RCI'   #  Vapor->Water
    
                Processes_signed_colors['PRR_WAU']='salmon'   #  Vapor->Water
                Processes_signed_names['PRR_WAU']='PRR_WAU'   #  Vapor->Water
    
                Processes_signed_colors['PRR_RCS']='cyan'   #  Vapor->Water
                Processes_signed_names['PRR_RCS']='PRR_RCS'   #  Vapor->Water
    
                Processes_signed_colors['PRR_RCG']='khaki'   #  Rain->Graupel
                Processes_signed_names['PRR_RCG']='PRR_RCG'   #  Rain->Graupel
    
                Processes_signed_colors['PRI_IHM']='darkorange'   #  Vapor->Water
                Processes_signed_names['PRI_IHM']='PRI_IHM'   #  Vapor->Water
    
                Processes_signed_colors['PRI_RFZ']='orange'   #  Vapor->Water
                Processes_signed_names['PRI_RFZ']='PRI_RFZ'   #  Vapor->Water
    
                Processes_signed_colors['PRI_RCI']='darkblue'   #  Vapor->Water
                Processes_signed_names['PRI_RCI']='PRI_RCI'   #  Vapor->Water
    
                Processes_signed_colors['PRI_IHA']='azure'   #  Vapor->Water
                Processes_signed_names['PRI_IHA']='PRI_IHA'   #  Vapor->Water
    
                Processes_signed_colors['PRS_SCI']='coral'   #  Vapor->Water
                Processes_signed_names['PRS_SCI']='PRS_SCI'   #  Vapor->Water
    
                Processes_signed_colors['PRS_RCS']='mediumseagreen'   #  Vapor->Water
                Processes_signed_names['PRS_RCS']='PRS_RCS'   #  Vapor->Water
    
                Processes_signed_colors['PRS_SCW']='mediumpurple'   #  Vapor->Water
                Processes_signed_names['PRS_SCW']='PRS_SCW'   #  Vapor->Water
    
                Processes_signed_colors['PRS_IHM']='skyblue'   #  Vapor->Water
                Processes_signed_names['PRS_IHM']='PRS_IHM'   #  Vapor->Water
    
                Processes_signed_colors['PRS_IDE']='beige'   #  Vapor->Water
                Processes_signed_names['PRS_IDE']='PRS_IDE'   #  Vapor->Water
    
                Processes_signed_colors['PRG_SCW']='violet'   #  Vapor->Water
                Processes_signed_names['PRG_SCW']='PRG_SCW'   #  Vapor->Water
    
                Processes_signed_colors['PRG_RFZ']='pink'   #  Vapor->Water
                Processes_signed_names['PRG_RFZ']='PRG_RFZ'   #  Vapor->Water
    
                Processes_signed_colors['PRG_GCW']='lightcyan'   #  Vapor->Water
                Processes_signed_names['PRG_GCW']='PRG_GCW'   #  Vapor->Water
    
                Processes_signed_colors['PRG_RCI']='lime'   #  Vapor->Water
                Processes_signed_names['PRG_RCI']='PRG_RCI'   #  Vapor->Water
    
                Processes_signed_colors['PRG_RCS']='peru'   #  Vapor->Water
                Processes_signed_names['PRG_RCS']='PRG_RCS'   #  Vapor->Water
    
                Processes_signed_colors['PRG_IHM']='violet'   #  Graupel->Ice
                Processes_signed_names['PRG_IHM']='PRG_IHM'   #  Graupel->Ice
    
                Processes_signed_colors['E_PRR_RCS']='pink'   #  Graupel->Ice
                Processes_signed_names['E_PRR_RCS']='E_PRR_RCS'   #  Graupel->Ice
    
                Processes_signed_colors['E_PRR_RCG']='yellow'   #  Graupel->Ice
                Processes_signed_names['E_PRR_RCG']='E_PRR_RCG'   #  Graupel->Ice
    
                Processes_signed_colors['E_PRI_IDE']='red'   #  Graupel->Ice
                Processes_signed_names['E_PRI_IDE']='E_PRI_IDE'   #  Graupel->Ice
    
                Processes_signed_colors['E_PRS_RCS']='green'   #  Graupel->Ice
                Processes_signed_names['E_PRS_RCS']='E_PRS_RCS'   #  Graupel->Ice
    
                Processes_signed_colors['E_PRG_GDE']='gold'   #  Graupel->Ice
                Processes_signed_names['E_PRG_GDE']='E_PRG_GDE'   #  Graupel->Ice
    
                Processes_signed_colors['PRI_WFI']='maroon'   #  Graupel->Ice
                Processes_signed_names['PRI_WFI']='PRI_WFI'   #  Graupel->Ice
                
                Processes_signed_colors['PRW_IMI']='maroon'   #  Graupel->Ice
                Processes_signed_names['PRW_IMI']='PRW_IMI'   #  Graupel->Ice


            
    return(Processes_signed_colors,Processes_signed_names)

color_condensation=   '#4b86c2'   #  bright blue
color_evaporation=    '#ffd633'   #  yellow
color_freezing=       '#002266'   #  dark blue
color_melting=        '#ff9a00'   #  orange
color_deposition=     '#34dabe'   #  cyan
color_sublimation=    '#8d5cbf'    #  purple
color_autoconversion= '#90db26'   #  green
color_ice=            '#D4F0FF'   #  icy blue




list_lumped_names_morrison=[]
list_lumped_processes_morrison=[]
lumped_colors_morrison={}

list_lumped_names_morrison.append('Condensation')
list_lumped_processes_morrison.append(['PCC','PCCN'])
lumped_colors_morrison['Condensation']=color_condensation

list_lumped_names_morrison.append('Evaporation')
list_lumped_processes_morrison.append(['EPCC','PRE'])
lumped_colors_morrison['Evaporation']=color_evaporation

list_lumped_names_morrison.append('Freezing')
list_lumped_processes_morrison.append(['MNUCCC','MNUCCR','PSACWS','PSACWI','PIACR','QMULTG','QMULTRG','QMULTS','QMULTR','PSACWG','PGRACS','PGSACW','QICF','QGRF','QNIRF'])
lumped_colors_morrison['Freezing']=color_freezing

list_lumped_names_morrison.append('Melting')
list_lumped_processes_morrison.append(['PSMLT','PGMLT','PRACG','PRACS','PRACGI','PSACR','PSACG','QIIM'])
lumped_colors_morrison['Melting']=color_melting

list_lumped_names_morrison.append('Autoconv./Accr.')
list_lumped_processes_morrison.append(['PRC','PRA'])
lumped_colors_morrison['Autoconv./Accr.']=color_autoconversion

list_lumped_names_morrison.append('Deposition')
list_lumped_processes_morrison.append(['PRD','PRDS','PRDG','MNUCCD'])
lumped_colors_morrison['Deposition']=color_deposition

list_lumped_names_morrison.append('Sublimation')
list_lumped_processes_morrison.append(['EPRD','EPRDS','EPRDG','EVPMS','EVPMG',])
lumped_colors_morrison['Sublimation']=color_sublimation

list_lumped_names_morrison.append('Ice processes')
list_lumped_processes_morrison.append(['PRAI','EPRCI','PRCI','PRACIS'])
lumped_colors_morrison['Ice processes']=color_ice

lumped_colors_morrison['Other']='grey'
                
                      
list_lumped_names_thompson=[]
list_lumped_processes_thompson=[]
lumped_colors_thompson={}

list_lumped_names_thompson.append('Condensation')
list_lumped_processes_thompson.append(['E_PRW_VCD'])
lumped_colors_thompson['Condensation']=color_condensation

list_lumped_names_thompson.append('Evaporation')
list_lumped_processes_thompson.append(['PRW_VCD','PRV_REV'])
lumped_colors_thompson['Evaporation']=color_evaporation

list_lumped_names_thompson.append('Freezing')
list_lumped_processes_thompson.append(['PRG_RFZ','PRI_WFZ','PRI_RFZ','PRG_RCG','PRG_GCW','PRG_SCW','PRS_SCW','PRS_RCS','PRG_RCS','PRR_RCI','PRI_WFI'])
lumped_colors_thompson['Freezing']=color_freezing


list_lumped_names_thompson.append('Melting')
list_lumped_processes_thompson.append(['PRR_SML','PRR_GML','PRR_RCG','PRR_RCS','PRR_RCS', 'PRW_IMI'])
lumped_colors_thompson['Melting']=color_melting

list_lumped_names_thompson.append('Autoconv./Accr.')
list_lumped_processes_thompson.append(['PRR_WAU','PRR_RCW'])
lumped_colors_thompson['Autoconv./Accr.']=color_autoconversion

list_lumped_names_thompson.append('Deposition')
list_lumped_processes_thompson.append(['PRG_GDE','PRS_SDE','PRI_IDE','PRS_IDE','PRI_INU','PRI_IHA'])
lumped_colors_thompson['Deposition']=color_deposition

list_lumped_names_thompson.append('Sublimation')
list_lumped_processes_thompson.append(['E_PRS_SDE','E_PRG_GDE','E_PRI_SDI'])
lumped_colors_thompson['Sublimation']=color_sublimation

list_lumped_names_thompson.append('Ice processes')
list_lumped_processes_thompson.append(['PRI_IHA','PRS_SCI','PRS_IAU','PRI_IHM','PRS_IHM','PRG_IHM'])
lumped_colors_thompson['Ice processes']=color_ice


lumped_colors_thompson['Other']='grey'




def lump_cubelist(cubelist_in,list_names_in, list_cubes_in,lumping='all',others=True):
    from iris.cube import CubeList
    if lumping=='all':
        list_names=list_names_in
        list_cubes=list_cubes_in

    if lumping=='latent':
        list_names=[]
        list_cubes=[]
        for i,name in enumerate(list_names_in):
            if name in ['Condensation','Evaporation','Freezing','Melting','Deposition','Sublimation']:
                list_names.append(list_names_in[i])
                list_cubes.append(list_cubes_in[i])

    if lumping=='mass':
        list_names=[]
        list_cubes=[]
        for i,name in enumerate(list_names_in):
            if name in ['Autoconv./Accr.','Ice processes']:
                list_names.append(list_names_in[i])
                list_cubes.append(list_cubes_in[i])


    
    cubelist_out=CubeList()
    list_cubes_other=[cube.name() for cube in cubelist_in]
    for i,name in enumerate(list_names):
        # Summ al cubes in list_cubes[i] and add them to output cubelist:
        #print(cubelist_in)
        #print(list_cubes[i])

        cubelist=cubelist_in.extract(list_cubes[i])
        if cubelist:
            cube=sum(cubelist)
            cube.rename(name)        
            cubelist_out.append(cube)
        #Remove these list_cubes from "Other"
        list_cubes_other=list(set(list_cubes_other)-set(list_cubes[i]))
    #Addd allremaining list_cubes and call them "other"
    if others:
        cubelist=cubelist_in.extract(list_cubes_other)
        if cubelist:
            cube=sum(cubelist)
            cube.rename('Other')
            cubelist_out.append(cube)

    return cubelist_out


def lump_processes(processes_in,microphysics=None,lumping='all',others=True):
    if (microphysics=='morrison'):
        processes_out=lump_cubelist(processes_in,list_lumped_names_morrison, list_lumped_processes_morrison,lumping=lumping,others=others)
    elif (microphysics=='thompson'):
        processes_out=lump_cubelist(processes_in,list_lumped_names_thompson, list_lumped_processes_thompson,lumping=lumping,others=others)       
    else:
        print('microphysics must be morrison or thompson')
    return(processes_out)



#def load_wrf_latenheating:
    
    
    