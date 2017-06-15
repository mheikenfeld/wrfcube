def load(filenames,variable,mode='auto',**kwargs):
    if mode=='auto':
        variable_list_file=variable_list(filenames)
        if variable in variable_list_file:
            variable_cube=loadwrfcube(filenames,variable,**kwargs)
        elif variable in variable_list_derive:
            variable_cube=derivewrfcube(filenames,variable,**kwargs)
        elif variable in variable_dict_pseudonym.keys:
            variable_load=variable_dict_pseudonym[variable]
            variable_cube=loadwrfcube(filenames,variable_load,**kwargs)

            
        
    elif mode=='file':
        variable_list_file=variable_list(filenames)
        if variable in variable_list_file:
            variable_cube=loadwrfcube(filenames,variable,**kwargs)
    elif mode=='derive':
        variable_cube=derivewrfcube(filenames,variable,**kwargs)
    elif mode=='pseudonym':
        variable_load=variable_dict_pseudonym[variable]
        variable_cube=loadwrfcube(filenames,variable_load,**kwargs)
    else:
       print("mode=",mode)
       raise SystemExit('unknown mode')

    return variable_cube




def loadwrfcubelist(filenames,variable_list,**kwargs):
    from iris.cube import CubeList
    cubelist_out=CubeList()  
    for variable in variable_list:
        cubelist_out.append(loadwrfcube(filenames,variable,**kwargs))
    return(cubelist_out)


def loadwrfcube(filenames,variable,**kwargs):
#    print(' in loadwrfcube: filenames= ',filenames)
#    print(' in loadwrfcube: variable= ',variable)
    if 'lazy' in kwargs:
        lazy=kwargs.pop('lazy')
    else:
        lazy=True

    if type(filenames) is list:
        variable_cube=loadwrfcube_mult(filenames,variable,**kwargs)
    elif type(filenames) is str:
        variable_cube=loadwrfcube_single(filenames,variable,**kwargs)
    else:
        print("filenames=",filenames)
        raise SystemExit('Type of input unknown: Must be str of list')
    
    # load data to get around bugs in lazy evaluation:
    if not lazy:
        variable_cube_data=variable_cube.data
    
    if 'add_coordinates' in kwargs:
        add_coordinates=kwargs['add_coordinates']
    else:
        add_coordinates=None
    if add_coordinates != None:
        add_aux_coordinates_multidim(filenames,variable_cube,**kwargs) 
    return variable_cube
    
def loadwrfcube_single(filenames,variable,constraint=None,add_coordinates=None):
    from iris import load_cube 
    variable_cube=load_cube(filenames,variable)
    
    variable_cube=addcoordinates(filenames, variable,variable_cube,add_coordinates=add_coordinates)        
    variable_cube=variable_cube.extract(constraint)
    return variable_cube
        
    
def loadwrfcube_mult(filenames,variable,constraint=None,add_coordinates=None):
    from iris.cube import CubeList
    #print(' in loadwrfcube_mult:',constraint)
    cube_list=[]
#    print(' in loadwrfcube_mult: filenames= ',filenames)
#    print(' in loadwrfcube_mult: variable= ',variable)
    for i in range(len(filenames)):
        cube_list.append(loadwrfcube_single(filenames[i],variable,add_coordinates=add_coordinates) )
    for member in cube_list:
        member.attributes={}
    variable_cubes=CubeList(cube_list)
    variable_cube=variable_cubes.concatenate_cube()
   # print(variable)
    #print(variable_cube)
    variable_cube=variable_cube.extract(constraint)
    return variable_cube

    
#def loadwrfcube_dimcoord(filenames,variable,**kwargs):
#    from iris import load_cube     
#    variable_cube=load_cube(filenames,variable)
#    for coord in variable_cube.coords():
#        variable_cube.remove_coord(coord.name())
#    variable_cube=add_dim_coordinates(filenames, variable,variable_cube,**kwargs)
#    #variable_cube=variable_cube.extract(constraint)
##    variable_cube_data=variable_cube.data
#    return variable_cube

    
#def loadwrfcube_nocoord(filenames,variable):
#    from iris import load_cube   
#    variable_cube=load_cube(filenames,variable)
#    for coord in variable_cube.coords():
#        variable_cube.remove_coord(coord.name())
##    variable_cube_data=variable_cube.data
#    return variable_cube

def derivewrfcubelist(filenames,variable_list,**kwargs):
    from iris.cube import CubeList
    cubelist_out=CubeList()  
    for variable in variable_list:
        cubelist_out.append(derivewrfcube(filenames,variable))
    return(cubelist_out)
    
#
#def derivewrfcube(filenames,variable,**kwargs):
#    if type(filenames) is list:
#        variable_cube=derivewrfcube_mult(filenames,variable,**kwargs)
#    elif type(filenames) is str:
#        variable_cube=derivewrfcube_single(filenames,variable,**kwargs)
#    else:
#        print('Type of input unknown: Must be str of list')
#    return variable_cube
#
#def derivewrfcube_mult(filenames,variable,**kwargs):
#    from iris.cube import CubeList
#    cube_list=[]
#    for i in range(len(filenames)):
#        cube_list.append(derivewrfcube_single(filenames[i],variable,**kwargs) )
#    for member in cube_list:
#        member.attributes={}
#    variable_cubes=CubeList(cube_list)
#    variable_cube=variable_cubes.concatenate_cube()
#    #variable_cube=variable_cube.extract(**kwargs.pop('constraint'))
#
#    return variable_cube

variable_dict_pseudonym={}
variable_dict_pseudonym['radar_relfectivity']='REFL10CM'




variable_list_derive=[
        'potential_temperature',
        'temperature','air_temperature',
        'density'
        'LWC',
        'IWC',
        'LWP',
        'IWP',
        'IWV',
        'airmass',
        'layer_height',
        'geopotential_height',
        'pressure',
        'relative_humidity',
        'w_at_T',
        'maximum_reflectivity' ,
        'surface_precipitation'
        ]


#def derivewrfcube_single(filenames,variable,**kwargs):
def derivewrfcube(filenames,variable,**kwargs):
    if variable == 'potential temperature':
        variable_cube=calculate_wrf_potential_temperature(filenames,**kwargs)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable in ['temperature','air_temperature']:
        variable_cube=calculate_wrf_temperature(filenames,**kwargs)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'density':
        variable_cube=calculate_wrf_density(filenames,**kwargs)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'LWC':    
        variable_cube=calculate_wrf_LWC(filenames,**kwargs)
        #variable_cube=addcoordinates(filenames, 'QCLOUD',variable_cube,add_coordinates)
    elif variable == 'IWC':  
        variable_cube=calculate_wrf_IWC(filenames,**kwargs)
        #variable_cube=addcoordinates(filenames, 'QICE',variable_cube,add_coordinates)    
    elif variable == 'LWP':    
        variable_cube=calculate_wrf_LWP(filenames,**kwargs)
        #variable_cube=addcoordinates(filenames, 'OLR',variable_cube,add_coordinates)
    elif variable == 'IWP':    
        variable_cube=calculate_wrf_IWP(filenames,**kwargs)
        #variable_cube=addcoordinates(filenames, 'OLR',variable_cube,add_coordinates)
    elif variable == 'IWV':    
        variable_cube=calculate_wrf_IWV(filenames,**kwargs)
        #variable_cube=addcoordinates(filenames, 'OLR',variable_cube,add_coordinates)
    elif variable == 'airmass':    
        variable_cube=calculate_wrf_airmass(filenames,**kwargs)
    elif variable == 'layer_height':    
        variable_cube=calculate_wrf_layerheight(filenames,**kwargs)
    elif variable == 'area':    
        variable_cube=calculate_wrf_area(filenames,**kwargs)        
    elif variable == 'geopotential_height':    
        variable_cube=calculate_wrf_geopotential_height(filenames,**kwargs)
        replace_cube=loadwrfcube(filenames,'T',**kwargs)
        variable_cube=replacecoordinates(variable_cube,replace_cube)  
    
    elif variable == 'geopotential_height_stag':    
        variable_cube=calculate_wrf_geopotential_height_stag(filenames,**kwargs)

#    elif variable == 'geopotential_height_xstag':    
#        variable_cube=calculate_wrf_geopotential_height_xstag(filenames,**kwargs)
#        replace_cube=loadwrfcube(filenames,'U',**kwargs)
#        variable_cube=replacecoordinates(variable_cube,replace_cube)  
#
#    elif variable == 'geopotential_height_ystag':    
#        variable_cube=calculate_wrf_geopotential_height_ystag(filenames,**kwargs)
#        replace_cube=loadwrfcube(filenames,'V',**kwargs)
#        variable_cube=replacecoordinates(variable_cube,replace_cube)  

    elif variable == 'pressure':    
        variable_cube=calculate_wrf_pressure(filenames,**kwargs)
        
    elif variable == 'geopotential':    
        variable_cube=calculate_wrf_geopotential(filenames,**kwargs)

    elif variable == 'pressure_xstag':    
        variable_cube=calculate_wrf_pressure(filenames,**kwargs)
        replace_cube=loadwrfcube(filenames,'U',**kwargs)
        variable_cube=replacecoordinates(variable_cube,replace_cube)  

    elif variable == 'pressure_ystag':    
        variable_cube=calculate_wrf_pressure(filenames,**kwargs)
        replace_cube=loadwrfcube(filenames,'V',**kwargs)
        variable_cube=replacecoordinates(variable_cube,replace_cube)  

    elif variable == 'relative_humidity':    
        variable_cube=calculate_wrf_relativehumidity(filenames,**kwargs)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'w_at_T':    
        variable_cube=calculate_wrf_w_at_T(filenames,**kwargs)
        replace_cube=loadwrfcube(filenames,'T',**kwargs)
        variable_cube=replacecoordinates(variable_cube,replace_cube)        
    elif variable == 'surface_precipitation':
        variable_cube=calculate_wrf_surface_precipitation(filenames,**kwargs)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'maximum reflectivity':    
        variable_cube=calculate_wrf_maximum_reflectivity(filenames,**kwargs)
    else:
        raise NameError(variable, 'is not a known variable') 
    return variable_cube
    
def calculate_wrf_surface_precipitation(filenames,constraint=None,add_coordinates=None):
    import numpy as np
    RAINNC= loadwrfcube(filenames, 'RAINNC',add_coordinates=add_coordinates)
    dt=(RAINNC.coords('time')[0].points[1]-RAINNC.coords('time')[0].points[0])*24
    rainnc_inst=np.concatenate((RAINNC.data[[1],:,:]-RAINNC.data[[0],:,:],RAINNC.data[1:,:,:]-RAINNC.data[0:-1:,:,:]),axis=0)/dt
    RAINNC_inst=RAINNC
    RAINNC_inst.data=rainnc_inst
    RAINNC_inst.rename('surface precipitation')
    RAINNC_inst.units= 'mm/h'
    RAINNC_inst=RAINNC_inst.extract(constraint)
    return RAINNC_inst

def variable_list(filenames):
    from netCDF4 import Dataset
    if type(filenames)==list:
        filenames=filenames[0]
    variable_list = list(Dataset(filenames).variables)
    return variable_list

def calculate_wrf_potential_temperature(filenames,**kwargs):
    from iris import coords
    T= loadwrfcube(filenames, 'T',**kwargs)
    T0 = coords.AuxCoord(300.0,long_name='reference_Temperature', units='K')
    theta=T+T0;
    theta.rename('potential temperature')
    return theta
    
def calculate_wrf_temperature(filenames,**kwargs):
    from iris import coords
    theta= derivewrfcube(filenames,'potential temperature',**kwargs)  
    p = derivewrfcube(filenames,'pressure',**kwargs)
    p0 =coords.AuxCoord(1000.0,long_name='reference_pressure', units='hPa')
    p0.convert_units(p.units)
    p1=p/p0
    exp=(287.05 / 1005.0)
    T=theta*(p1**exp) #work around iris issue here by loading one of the cubes into numpy array..
    T.rename('air temperature')
    return T
    
def calculate_wrf_relativehumidity(filenames,**kwargs):
    from iris import cube
    QVAPOR=loadwrfcube(filenames, 'QVAPOR',**kwargs)
    T=derivewrfcube(filenames,'temperature',**kwargs)
    p=derivewrfcube(filenames,'pressure',**kwargs)
    p.convert_units('Pa')
    p=p.data
    rh=calculate_RH(QVAPOR.data,T.data,p.data)
    RH=cube.Cube(rh, units='percent',long_name='realtive humidity')
    return RH   
    
def calculate_RH(QVAPOR,T,p):
    from numpy import exp,maximum,minimum
    ES=1e2*6.1094*exp(17.625*(T-273.15)/(T-273.15+243.04))
    QVS = 0.622*ES/ (p- (1.0-0.622)*ES)
    RH = 100*maximum(minimum(QVAPOR/QVS,1.0),0.0)
    return RH   
    
def calculate_wrf_LWC(filenames,**kwargs):
    microphysics_scheme=kwargs.pop(microphysics_scheme)
    #QCLOUD=loadwrfcube(filenames, 'QCLOUD',**kwargs)
    #QRAIN=loadwrfcube(filenames, 'QRAIN',**kwargs)
    #LWC=QCLOUD+QRAIN    
    list_variables=['QCLOUD','QRAIN']
    LWC=load_sum(filename,list_variables,**kwargs)
    LWC.rename('liquid water content')
    #LWC.rename('mass_concentration_of_liquid_water_in_air')
    return LWC   
#    
def calculate_wrf_IWC(filenames,**kwargs):
    microphysics_scheme=kwargs.pop(microphysics_scheme)
    #QICE=loadwrfcube(filenames, 'QICE',**kwargs)
    #QSNOW=loadwrfcube(filenames, 'QSNOW',**kwargs)
    #QGRAUP=loadwrfcube(filenames, 'QGRAUP',**kwargs)
    #IWC=QICE+QSNOW+QGRAUP
    if microphysics_scheme in [None,morrison,thompson]:
        list_variables=['QICE','QSNOW','QGRAUP']
    elif microphysics_scheme in ["SBM_full"]:
        list_variables=['QICE','QSNOW','QGRAUP']
    elif microphysics_scheme in ["SBM_full"]:
        list_variables=['QICEC','QICED','QICEP','QSNOW','QGRAUP','QHAIL']
    IWC=load_sum(filename,list_variables,**kwargs)
    IWC.rename('ice water content')
    #IWC.rename('mass_concentration_of_ice_water_in_air')

    return IWC
    
def calculate_wrf_airmass(filenames,**kwargs):
    rho=derivewrfcube(filenames,'density',**kwargs)
    layer_height=derivewrfcube(filenames,'layer_height',**kwargs)
    Airmass=rho*layer_height
    Airmass.rename('mass of air')
    return Airmass
    
    
def calculate_wrf_layerheight(filenames,**kwargs):
    from iris import Constraint
    zH=derivewrfcube(filenames,'geopotential_height_stag',**kwargs)
    bottom_top_stag=zH.coord('bottom_top_stag').points
    layer_height = (zH.extract(Constraint(bottom_top_stag=bottom_top_stag[1:]))-zH.extract(Constraint(bottom_top_stag=bottom_top_stag[:-1])).data) 
    layer_height.rename('layer_height')
    replace_cube=loadwrfcube(filenames,'T',**kwargs)
    layer_height=replacecoordinates(layer_height,replace_cube)
    return layer_height
    
def calculate_wrf_LWP(filenames,**kwargs):
    from iris.analysis import SUM
    LWC=derivewrfcube(filenames,'LWC',**kwargs)
    Airmass=derivewrfcube(filenames,'airmass',**kwargs)
    LWP=(LWC*Airmass).collapsed(('bottom_top'),SUM)
    LWP.rename('liquid water path')
    #LWP.rename('atmosphere_mass_content_of_cloud_liquid_water')
    return LWP   
#    
def calculate_wrf_IWP(filenames,**kwargs):    
    from iris.analysis import SUM
    IWC=derivewrfcube(filenames,'IWC',**kwargs)
    Airmass=derivewrfcube(filenames,'airmass',**kwargs)
    IWP=(IWC*Airmass).collapsed(('bottom_top'),SUM)
    IWP.rename('ice water path')
    #IWP.rename('atmosphere_mass_content_of_cloud_ice_water')
    return IWP
    
def calculate_wrf_IWV(filenames,**kwargs):    
    from iris.analysis import SUM
    QVAPOR=loadwrfcube(filenames,'QVAPOR',**kwargs)
    Airmass=derivewrfcube(filenames,'airmass',**kwargs)
    IWV=(QVAPOR*Airmass).collapsed(('bottom_top'),SUM)
    IWV.rename('integrated water vapour')
    #IWV.rename('atmosphere_mass_content_of_water_vapor')
    return IWV

def integrate_cube(variable,Airmass_or_dz):
    from iris.analysis import SUM
    name=variable.name()
    variable_integrated=(variable*Airmass_or_dz)
    variable_integrated.remove_coord('geopotential_height')
    variable_integrated=variable_integrated.collapsed(('bottom_top'),SUM)
    variable_integrated.rename(name)
    return variable_integrated
    
def calculate_wrf_LWP_fromcubes(LWC,Airmass):
    from iris.analysis import SUM
    LW=(LWC*Airmass)
    LW.remove_coord('geopotential_height')
    LWP=LW.collapsed(('bottom_top'),SUM)
    LWP.rename('liquid water path')
    #LWP.rename('atmosphere_mass_content_of_cloud_liquid_water')
    return LWP
    
#    
def calculate_wrf_IWP_fromcubes(IWC, Airmass):
    from iris.analysis import SUM
    IW=(IWC*Airmass)
    IW.remove_coord('geopotential_height')
    IWP=IW.collapsed(('bottom_top'),SUM)
    IWP.rename('ice water path')
    return IWP

def calculate_wrf_IWV_fromcubes(QVAPOR,Airmass):
    from iris.analysis import SUM
    VAPOR=(QVAPOR*Airmass)
    VAPOR.remove_coord('geopotential_height')
    IWV=VAPOR.collapsed(('bottom_top'),SUM)
    IWV.rename('integrated water vapor')
    return IWV


    
def calculate_wrf_maximum_reflectivity(filenames,**kwargs):
    from iris.analysis import MAX
    REFL_10CM=loadwrfcube(filenames,'REFL_10CM',**kwargs)
    MAX_REFL_10CM=REFL_10CM.collapsed('bottom_top', MAX)
    MAX_REFL_10CM.rename('maximum reflectivity')
    return MAX_REFL_10CM
    
def calculate_wrf_w_at_T(filenames,**kwargs):
    from iris import cube,Constraint
    w=loadwrfcube(filenames, 'W',**kwargs)
    constraint_1=Constraint(bottom_top_stag=lambda cell:  cell > w.coord('bottom_top_stag').points[0])
    constraint_2=Constraint(bottom_top_stag=lambda cell:  cell < w.coord('bottom_top_stag').points[-1])
    w_at_T = cube.Cube(0.5*(w.extract(constraint_1).data+w.extract(constraint_2).data),var_name='w',long_name='vertical velocity on T grid', units='m/s')

    return w_at_T

def calculate_wrf_density(filenames,**kwargs):
    from iris import coords

    if ('ALT' in variable_list(filenames)):
        alt=loadwrfcube(filenames,'ALT',**kwargs)
        rho=alt**(-1)
    else: 
       R=coords.AuxCoord(287.058,long_name='Specific gas constant for air',units='Joule kg^-1 K^-1')
       p=derivewrfcube(filenames,'pressure',**kwargs)
       T=derivewrfcube(filenames,'temperature',**kwargs)
       rho=p*((R*T)**-1)
       rho.rename('air_density')
    return rho
#    
def calculate_wrf_pressure(filenames,**kwargs):
    P= loadwrfcube(filenames, 'P',**kwargs)
    PB= loadwrfcube(filenames, 'PB',**kwargs)
    p=P + PB 
    p.rename('pressure')
    return p
    

def calculate_wrf_pressure_stag(filenames,**kwargs):
    from iris import Constraint
    p= derivewrfcube(filenames, 'pressure',**kwargs)
    bottom_top=p.coord('bottom_top').points
    p_stag = 0.5*(p.extract(Constraint(bottom_top=bottom_top[:-1]))+p.extract(Constraint(bottom_top=bottom_top[1:])).data)    
    return p_stag

    #    
def calculate_wrf_pressure_xstag(filenames,**kwargs):
    from iris import Constraint
    p=derivewrfcube(filenames,'pressure',**kwargs)
    west_east=p.coord('west_east').points
    p_xstag = 0.5*(p.extract(Constraint(west_east=west_east[:-1]))+p.extract(Constraint(west_east=west_east[1:])).data)    
    p_xstag.rename('pressure')
    return p

#    
def calculate_wrf_pressure_ystag(filenames,**kwargs):
    from iris import Constraint
    p=derivewrfcube(filenames,'pressure',**kwargs)
    south_north=p.coord('south_north').points
    p_ystag = 0.5*(p.extract(Constraint(south_north=south_north[:-1]))+p.extract(Constraint(south_north=south_north[1:])).data)
    p_ystag.rename('pressure')
    return p_ystag

    
def calculate_wrf_geopotential(filenames,**kwargs):
    PH= loadwrfcube(filenames, 'PH',**kwargs)
    PHB= loadwrfcube(filenames,'PHB',**kwargs)
    pH=PH + PHB
    pH.rename('geopotential')
    return pH
    
def unstagger(cube_in,coord,filenames,**kwargs):
    cube_out=cube_interp_reduceby1(cube_in,coord)
    replace_cube=loadwrfcube(filenames,'T',**kwargs)
    cube_out=replacecoordinates(cube_out,replace_cube)


def array_interp_extendby1(array,dim):
    import numpy as np
    idx1=[slice(None)] * (array.ndim)
    idx2=[slice(None)] * (array.ndim)
    idx_start=[slice(None)] * (array.ndim)
    idx_end=[slice(None)] * (array.ndim)
    idx1[dim]=slice(1,None)
    idx2[dim]=slice(0,-1)
    idx_start[dim]=slice(1)
    idx_end[dim]=slice(-2,-1)
    array_out=np.concatenate((array[idx_start],0.5*(array[idx1]+array[idx2]),array[idx_end]),axis=dim)
    return array_out
    
def array_interp_reduceby1(array,dim):
    idx = [slice(None)] * (array.ndim)  
    idx1=idx
    idx2=idx
    idx1[dim]=slice(1,None)
    idx2[dim]=slice(None,-1)

    array_out=0.5*(array[idx1]+array[idx2])
    return array_out


def cube_interp_extendby1(cube_in,coord):
    import numpy as np
    dim=cube_in.coord_dims(coord)[0]
    cube_data=cube_in.data
    ndim=cube_in.ndim
    idx1=[slice(None)] * (ndim)
    idx2=[slice(None)] * (ndim)
    idx_start=[slice(None)] * (ndim)
    idx_end=[slice(None)] * (ndim)
    idx1[dim]=slice(1,None)
    idx2[dim]=slice(0,-1)
    idx_start[dim]=slice(1)
    idx_end[dim]=slice(-2,-1)
    cube_data[idx_start]
    cube_data[idx1]
    cube_data[idx2]
    cube_data[idx_end]    
    array_out=np.concatenate((cube_data[idx_start],0.5*(cube_data[idx1]+cube_data[idx2]),cube_data[idx_end]),axis=dim)

#    idx1[dim]=slice(1)
#    array_out=np.concatenate((cube_data[idx1],cube_data),axis=dim)
    return array_out
    
def cube_interp_reduceby1(cube_in,coord):
    dim=cube_in.coord_dims(coord)
    ndim=cube_in.ndim
    idx1=[slice(None)] * (ndim)
    idx2=[slice(None)] * (ndim)
    idx1[dim]=slice(1,None)
    idx2[dim]=slice(None,-1)

    cube_out=0.5*(cube_in[idx1]+cube_in[idx2].data)
    return cube_out

def load_sum(filename,list_variables,**kwargs):
    cube_out=load(filename,list_variables[0],**kwargs)
    for variable in list_variables[1:]:
        cube_out=cube_out+load(filename,variable,**kwargs)
    return cube_out
    


def calculate_wrf_geopotential_height_stag(filenames,**kwargs):
    from iris import coords
    pH=derivewrfcube(filenames,'geopotential',**kwargs)
    g = coords.AuxCoord(9.81,long_name='acceleration', units='m s^-2')
    zH=pH/g
    zH.rename('geopotential_height')
    return zH

def calculate_wrf_geopotential_height(filenames,**kwargs):
    from iris import Constraint
    zH=derivewrfcube(filenames,'geopotential_height_stag',**kwargs)
    bottom_top_stag=zH.coord('bottom_top_stag').points
    z = 0.5*(zH.extract(Constraint(bottom_top_stag=bottom_top_stag[:-1]))+zH.extract(Constraint(bottom_top_stag=bottom_top_stag[1:])).data) 
    z.rename('geopotential_height')
    #z = 0.5*(zH[:,:-1,:,:]+zH.data[:,1:,:,:])
    return z
    
def calculate_wrf_geopotential_height_ystag(filenames,**kwargs):
    #from iris import Constraint
    #z=derivewrfcube(filenames,'geopotential_height',**kwargs)
    z=calculate_wrf_geopotential_height(filenames,**kwargs)
    #south_north=z.coord('south_north').points
    #z_ystag = 0.5*(z.extract(Constraint(south_north=south_north[:-1]))+z.extract(Constraint(south_north=south_north[1:])).data)
    dim=z.coord_dims('south_north')[0]
    z_ystag = cube_interp_extendby1(z,'south_north')
    return z_ystag
    
def calculate_wrf_geopotential_height_xstag(filenames,**kwargs):
    #from iris import Constraint
   # z=derivewrfcube(filenames,'geopotential_height',**kwargs)
    z=calculate_wrf_geopotential_height(filenames,**kwargs)

    #west_east=z.coord('west_east').points
    #z_xstag = 0.5*(z.extract(Constraint(west_east=west_east[:-1]))+z.extract(Constraint(west_east=west_east[1:])).data)
    dim=z.coord_dims('west_east')[0]
    z_xstag = cube_interp_extendby1(z,'west_east')
    return z_xstag

def remove_all_coordinates(variable_cube):
    for coordinate in variable_cube.coords():
        variable_cube.remove_coord(coordinate.name())
    return variable_cube    

def replacecoordinates(variable_cube,replace_cube):        
    variable_cube_out=replace_cube
    variable_cube_out.data=variable_cube.data
    variable_cube_out.rename(variable_cube.name())
    variable_cube_out.units=variable_cube.units
    variable_cube_out.attributes={}
    return variable_cube_out

    
def addcoordinates(filenames, variable,variable_cube,**kwargs):
    if 'add_coordinates' in kwargs:
        add_coordinates=kwargs['add_coordinates']
    else:
        add_coordinates=None
    if add_coordinates==None:
        variable_cube=add_dim_coordinates(filenames, variable,variable_cube,**kwargs)
    else:
        variable_cube=add_dim_coordinates(filenames, variable,variable_cube,**kwargs)
        variable_cube=add_aux_coordinates_1dim(filenames, variable,variable_cube,**kwargs)
    return variable_cube

def add_dim_coordinates(filenames, variable,variable_cube,add_coordinates=None):
    from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
    from iris import load_cube
    variable_cube_dim= load_cube(filenames, variable)

    attributes=variable_cube_dim.attributes
    nc_id=Dataset(filenames)
    nc_variable=nc_id.variables[variable]
    variable_dimensions=nc_variable.dimensions
    [str(line) for line in variable_dimensions]
    DX=attributes['DX']        
    DY=attributes['DY']
    WEST_EAST_PATCH_END_UNSTAG=attributes['WEST-EAST_PATCH_END_UNSTAG']
    SOUTH_NORTH_PATCH_END_UNSTAG=attributes['SOUTH-NORTH_PATCH_END_UNSTAG']
    BOTTOM_TOP_PATCH_END_UNSTAG=attributes['BOTTOM-TOP_PATCH_END_UNSTAG']
    WEST_EAST_PATCH_END_STAG=attributes['WEST-EAST_PATCH_END_STAG']
    SOUTH_NORTH_PATCH_END_STAG=attributes['SOUTH-NORTH_PATCH_END_STAG']
    BOTTOM_TOP_PATCH_END_STAG=attributes['BOTTOM-TOP_PATCH_END_STAG']
    for dim in range(len(variable_dimensions)):
        if (variable_dimensions[dim]=='Time'):
           time=make_time_coord(filenames)
           variable_cube.add_dim_coord(time,dim)
        elif (variable_dimensions[dim]=='west_east'):
            west_east=make_westeast_coord(DX,WEST_EAST_PATCH_END_UNSTAG)
            variable_cube.add_dim_coord(west_east,dim)
        elif (variable_dimensions[dim]=='south_north'):
           south_north=make_southnorth_coord(DY, SOUTH_NORTH_PATCH_END_UNSTAG)
           variable_cube.add_dim_coord(south_north,dim)
        elif (variable_dimensions[dim]=='bottom_top'):
           bottom_top=make_bottom_top_coordinate(BOTTOM_TOP_PATCH_END_UNSTAG)   
           variable_cube.add_dim_coord(bottom_top,dim)
        elif variable_dimensions[dim]=='west_east_stag':
           west_east_stag=make_westeast_stag_coord(DX,WEST_EAST_PATCH_END_STAG)
           variable_cube.add_dim_coord(west_east_stag,dim)
        elif variable_dimensions[dim]=='south_north_stag':
           south_north_stag=make_southnorth_stag_coord(DY, SOUTH_NORTH_PATCH_END_STAG)
           variable_cube.add_dim_coord(south_north_stag,dim)
        elif variable_dimensions[dim]=='bottom_top_stag':
           bottom_top_stag=make_bottom_top_stag_coordinate(BOTTOM_TOP_PATCH_END_STAG)   
           variable_cube.add_dim_coord(bottom_top_stag,dim)
    return variable_cube        
    
def add_aux_coordinates_1dim(filenames, variable,variable_cube,add_coordinates=None):
    from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
    from iris import load_cube
    variable_cube_dim= load_cube(filenames, variable)
    attributes=variable_cube_dim.attributes
    nc_id=Dataset(filenames)
    
    nc_variable=nc_id.variables[variable]
    variable_dimensions=nc_variable.dimensions
    [str(line) for line in variable_dimensions]
    DX=attributes['DX']        
    DY=attributes['DY']
    WEST_EAST_PATCH_END_UNSTAG=attributes['WEST-EAST_PATCH_END_UNSTAG']
    SOUTH_NORTH_PATCH_END_UNSTAG=attributes['SOUTH-NORTH_PATCH_END_UNSTAG']
    #BOTTOM_TOP_PATCH_END_UNSTAG=attributes['BOTTOM-TOP_PATCH_END_UNSTAG']
    WEST_EAST_PATCH_END_STAG=attributes['WEST-EAST_PATCH_END_STAG']
    SOUTH_NORTH_PATCH_END_STAG=attributes['SOUTH-NORTH_PATCH_END_STAG']
    #BOTTOM_TOP_PATCH_END_STAG=attributes['BOTTOM-TOP_PATCH_END_STAG']
    coord_system=make_coord_system(attributes)
    coords=variable_cube.coords()
    if type(add_coordinates)!=list:
        add_coordinates1=add_coordinates
        add_coordinates=[]
        add_coordinates.append(add_coordinates1)
    for coordinate in add_coordinates:
        if coordinate=='xy':
            for dim in range(len(coords)):            
                if (coords[dim].name()=='west_east'):
                    x_coord=make_x_coord(DX,WEST_EAST_PATCH_END_UNSTAG,coord_system=coord_system)
                    variable_cube.add_aux_coord(x_coord,dim)
                elif (coords[dim].name()=='south_north'):
                    y_coord=make_y_coord(DY, SOUTH_NORTH_PATCH_END_UNSTAG,coord_system=coord_system)
                    variable_cube.add_aux_coord(y_coord,dim)
                elif (coords[dim].name()=='west_east_stag'):
                    x_stag_coord=make_x_stag_coord(DX,WEST_EAST_PATCH_END_STAG,coord_system=coord_system)
                    variable_cube.add_aux_coord(x_stag_coord,dim)
                elif coords[dim].name()=='south_north_stag':
                    y_stag_coord=make_y_stag_coord(DY, SOUTH_NORTH_PATCH_END_STAG,coord_system=coord_system)
                    variable_cube.add_aux_coord(y_stag_coord,dim)
    return variable_cube               
        
    
    
    
def add_aux_coordinates_multidim(filenames,variable_cube,**kwargs):
    import sys
    coords=variable_cube.coords()        
    add_coordinates=kwargs.pop('add_coordinates')
    if type(add_coordinates)!=list:
        add_coordinates1=add_coordinates
        add_coordinates=[]
        add_coordinates.append(add_coordinates1)
    for coordinate in add_coordinates:
     #   print(coordinate)

        if coordinate=='z':                
         #   print('coordinate is z')

            if (coords[0].name()=='time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east'): 
                z_coord=make_z_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(z_coord,(0,1,2,3))
            elif (coords[0].name()=='bottom_top' and coords[1].name()=='south_north' and coords[2].name()=='west_east'):
                z_coord=make_z_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(z_coord,(0,1,2))
            elif (coords[0].name()=='time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east_stag'):
                z_coord=make_z_xstag_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(z_coord,(0,1,2,3))
            elif (coords[0].name()=='bottom_top' and coords[1].name()=='south_north' and coords[2].name()=='west_east_stag'):
                z_coord=make_z_xstag_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(z_coord,(0,1,2))
            elif (coords[0].name()=='time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north_stag' and coords[3].name()=='west_east'):
                z_coord=make_z_ystag_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(z_coord,(0,1,2,3))
            elif (coords[0].name()=='bottom_top' and coords[1].name()=='south_north_stag' and coords[2].name()=='west_east'):
                z_coord=make_z_ystag_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(z_coord,(0,1,2))
            elif (coords[0].name()=='time' and coords[1].name()=='bottom_top_stag' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                z_stag_coord=make_z_stag_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(z_stag_coord,(0,1,2,3))
            elif (coords[0].name()=='bottom_top_stag' and coords[1].name()=='south_north' and coords[2].name()=='west_east'):
                z_stag_coord=make_z_stag_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(z_stag_coord,(0,1,2))
            else:
                print(coords)
                raise SystemExit("no z coordinates added")
                
        if coordinate=='pressure' :      
            if (coords[0].name()=='time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                p_coord=make_p_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            elif (coords[0].name()=='time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east_stag'):
                p_coord=make_p_xstag_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(p_coord)
            elif (coords[0].name()=='time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north_stag' and coords[3].name()=='west_east'):
                p_coord=make_p_ystag_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            elif (coords[0].name()=='time' and coords[1].name()=='bottom_top_stag' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                p_coord=make_p_stag_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            else:
                print(coords)
                raise SystemExit("p coordinates added")
                
        if (coordinate=='zp' or coordinate=='pz'):    
            if (coords[0].name()=='time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                z_coord=make_z_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(z_coord,(0,1,2,3))
                p_coord=make_p_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            elif (coords[0].name()=='time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east_stag'):
                z_coord=make_z_xstag_coordinate(filenames,**kwargs)     
                variable_cube.add_aux_coord(z_coord,(0,1,2,3))
                p_coord=make_p_xstag_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            elif (coords[0].name()=='time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north_stag' and coords[3].name()=='west_east'):
                z_coord=make_z_ystag_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(z_coord,(0,1,2,3))
                p_coord=make_p_ystag_coordinate(filenames,**kwargs)
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            elif (coords[0].name()=='time' and coords[1].name()=='bottom_top_stag' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                z_stag_coord=make_z_stag_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(z_stag_coord,(0,1,2,3))
                p_coord=make_p_stag_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            else:
                print(coords)
                raise SystemExit("no z and p coordinates added")
                
        if coordinate=='latlon':    
            if (coords[0].name()=='time' and (coords[1].name()=='bottom_top' or 'bottom_top_stag') and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                lat_coord=make_lat_coordinate(filenames,**kwargs)
                lon_coord=make_lon_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(lat_coord,(0,2,3))
                variable_cube.add_aux_coord(lon_coord,(0,2,3))            
            elif (coords[0].name()=='time' and (coords[1].name()=='bottom_top' or 'bottom_top_stag') and coords[2].name()=='south_north' and coords[3].name()=='west_east_stag'):
                lat_coord=make_lat_xstag_coordinate(filenames,**kwargs)
                lon_coord=make_lon_xstag_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(lat_coord,(0,2,3))
                variable_cube.add_aux_coord(lon_coord,(0,2,3))
            elif (coords[0].name()=='time' and (coords[1].name()=='bottom_top' or 'bottom_top_stag') and coords[2].name()=='south_north_stag' and coords[3].name()=='west_east'):
                lat_coord=make_lat_ystag_coordinate(filenames,**kwargs)
                lon_coord=make_lon_ystag_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(lat_coord,(0,2,3))
                variable_cube.add_aux_coord(lon_coord,(0,2,3))
            elif (coords[0].name()=='time'  and coords[1].name()=='south_north' and coords[2].name()=='west_east'):
                lat_coord=make_lat_coordinate(filenames,**kwargs)
                lon_coord=make_lon_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(lat_coord,(0,1,2))
                variable_cube.add_aux_coord(lon_coord,(0,1,2))            
            elif (coords[0].name()=='time'  and coords[1].name()=='south_north' and coords[2].name()=='west_east_stag'):
                lat_coord=make_lat_xstag_coordinate(filenames,**kwargs)
                lon_coord=make_lon_xstag_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(lat_coord,(0,1,2))
                variable_cube.add_aux_coord(lon_coord,(0,1,2))
            elif (coords[0].name()=='time' and coords[1].name()=='south_north_stag' and coords[2].name()=='west_east'):
                lat_coord=make_lat_ystag_coordinate(filenames,**kwargs)
                lon_coord=make_lon_ystag_coordinate(filenames,**kwargs)   
                variable_cube.add_aux_coord(lat_coord,(0,1,2))
                variable_cube.add_aux_coord(lon_coord,(0,1,2))
            else:
                print(coords)
                raise SystemExit("no lat/lon coordinates added")
    return variable_cube

    
def make_time_coord(filenames):
    from iris import load_cube,coords
    from datetime import datetime,timedelta
    from numpy import empty
    Times= load_cube(filenames, 'Times')
    filetimes = Times.data   
    filetimelist = []   # Will contain list of times in seconds since model start time in file.
    timeobjlist = []    # Will contain list of corresponding datetime objects
    for i, filetime in enumerate(filetimes):
        timeobj = datetime.strptime(filetime.tostring().decode('UTF-8'), "%Y-%m-%d_%H:%M:%S")
        if i == 0:
            timeobj0 = timeobj
        time_dt = timeobj-timeobj0 # timedelta object representing difference in time from start time
        filetimelist.append(time_dt.seconds)
        timeobjlist.append(timeobj)
    time_days=empty(len(timeobjlist))
    #Include a different base_date for dates close to 0001-01-01 (idealised simulations)
    if timeobjlist[0]<datetime(100,1,1):
        base_date=datetime(1,1,1)
    else:
        base_date=datetime(1970,1,1)
    time_units='days since '+ base_date.strftime('%Y-%m-%d')
    for i in range(len(timeobjlist)):
        time_days[i]=(timeobjlist[i] - base_date).total_seconds() / timedelta(1).total_seconds()
    time_coord=coords.DimCoord(time_days, standard_name='time', long_name='time', var_name='time', units=time_units, bounds=None, attributes=None, coord_system=None, circular=False)
    return time_coord              

def make_westeast_coord(DX,WEST_EAST_PATCH_END_UNSTAG):
    from iris import coords
    from numpy import arange
    WEST_EAST=arange(0,WEST_EAST_PATCH_END_UNSTAG)
    west_east=coords.DimCoord(WEST_EAST, standard_name=None, long_name='west_east', var_name='west_east', units='1', bounds=None, attributes=None, coord_system=None, circular=False)
    return west_east


def make_westeast_stag_coord(DX,WEST_EAST_PATCH_END_STAG):
    from iris import coords
    from numpy import arange
    WEST_EAST_U=arange(0,WEST_EAST_PATCH_END_STAG)
    west_east_stag=coords.DimCoord(WEST_EAST_U, standard_name=None, long_name='west_east_stag', var_name='west_east_stag', units='1', bounds=None, attributes=None, coord_system=None, circular=False)
    return west_east_stag

def make_southnorth_coord(DY,SOUTH_NORTH_PATCH_END_UNSTAG):
    from iris import coords
    from numpy import arange    #DY=attributes['DY']
    #SOUTH_NORTH_PATCH_END_UNSTAG=attributes['SOUTH_NORTH_PATCH_END_UNSTAG']
    SOUTH_NORTH=arange(0,SOUTH_NORTH_PATCH_END_UNSTAG)
    south_north=coords.DimCoord(SOUTH_NORTH, standard_name=None, long_name='south_north', var_name='south_north', units='1', bounds=None, attributes=None, coord_system=None, circular=False)
    return south_north
    
def make_southnorth_stag_coord(DY,SOUTH_NORTH_PATCH_END_STAG):
    from iris import coords
    from numpy import arange
    SOUTH_NORTH_V=arange(0,SOUTH_NORTH_PATCH_END_STAG)
    south_north_stag=coords.DimCoord(SOUTH_NORTH_V, standard_name=None, long_name='south_north_stag', var_name='south_north_stag', units='1', bounds=None, attributes=None, coord_system=None, circular=False)
    return south_north_stag

def make_bottom_top_coordinate(BOTTOM_TOP_PATCH_END_UNSTAG):
    from iris import coords
    from numpy import arange
    BOTTOM_TOP=arange(0,BOTTOM_TOP_PATCH_END_UNSTAG)
    bottom_top=coords.DimCoord(BOTTOM_TOP, standard_name=None, long_name='bottom_top', var_name='bottom_top', units='1', bounds=None, attributes=None, coord_system=None, circular=False)
    return bottom_top
    
def make_bottom_top_stag_coordinate(BOTTOM_TOP_PATCH_END_STAG):  
    from iris import coords
    from numpy import arange
    BOTTOM_TOP_W=arange(0,BOTTOM_TOP_PATCH_END_STAG)
    bottom_top_stag=coords.DimCoord(BOTTOM_TOP_W, standard_name=None, long_name='bottom_top_stag', var_name='bottom_top_stag', units='1', bounds=None, attributes=None, coord_system=None, circular=False)
    return bottom_top_stag

def make_coord_system(attributes):
    from iris import coord_systems
#    :CEN_LAT = -3.212929f ;
#		:CEN_LON = -60.59799f ;
#		:TRUELAT1 = 0.f ;
#		:TRUELAT2 = -5.f ;
#		:MOAD_CEN_LAT = -3.212929f ;
#		:STAND_LON = -60.f ;
#		:POLE_LAT = 90.f ;
#		:POLE_LON = 0.f ;
#		:GMT = 0.f ;
#		:JULYR = 2014 ;
#		:JULDAY = 244 ;
#		:MAP_PROJ = 1 ;
#		:MAP_PROJ_CHAR = "Lambert Conformal" ;
    MAP_PROJ_CHAR=attributes['MAP_PROJ_CHAR']    
    MAP_PROJ=attributes['MAP_PROJ']
    
    # cartesian coordinate system (idealized simulations):
    if (MAP_PROJ_CHAR=='Cartesian' and MAP_PROJ==0):
        coord_system=None
    
    # lambert Conformal system (idealized simulations):
    if (MAP_PROJ_CHAR=='Lambert Conformal' and MAP_PROJ==1):
        CEN_LON = attributes['CEN_LON']
        TRUELAT1 = attributes['TRUELAT1']
        TRUELAT2 = attributes['TRUELAT2']
        MOAD_CEN_LAT = attributes['MOAD_CEN_LAT']
        STAND_LON = attributes['STAND_LON']
        POLE_LAT = attributes['POLE_LAT']
        POLE_LON = attributes['POLE_LON']
        coord_system=coord_systems.LambertConformal(central_lat=MOAD_CEN_LAT, central_lon=CEN_LON, false_easting=0.0, false_northing=0.0, secant_latitudes=(TRUELAT1, TRUELAT2))
    return coord_system


def make_x_coord(DX,WEST_EAST_PATCH_END_UNSTAG,coord_system):
    from iris import coords
    from numpy import arange
    X=DX*(arange(0,WEST_EAST_PATCH_END_UNSTAG)-0.5)
    x_coord=coords.AuxCoord(X, standard_name='projection_x_coordinate', long_name='x', var_name='x', units='m', bounds=None, attributes=None, coord_system=coord_system)
    #x_coord.add_dim_coord(west_east,0)
    return x_coord

def make_x_stag_coord(DX,WEST_EAST_PATCH_END_STAG,coord_system=None):
    from iris import coords
    from numpy import arange
    X_U=DX*(arange(0,WEST_EAST_PATCH_END_STAG)-1)
    x_stag_coord=coords.AuxCoord(X_U, standard_name='projection_x_coordinate', long_name='x', var_name='x', units='m', bounds=None, attributes=None, coord_system=coord_system)
    #x_stag_coord.add_dim_coord(west_east_stag,0)
    return x_stag_coord

def make_y_coord(DY,SOUTH_NORTH_PATCH_END_UNSTAG,coord_system=None):
    from iris import coords
    from numpy import arange    #DY=attributes['DY']
    Y=DY*(arange(0,SOUTH_NORTH_PATCH_END_UNSTAG)-0.5)
    y_coord=coords.AuxCoord(Y, standard_name='projection_y_coordinate', long_name='y', var_name='y', units='m', bounds=None, attributes=None, coord_system=coord_system)
    #y_coord.add_dim_coord(south_north,0)
    return y_coord
    
def make_y_stag_coord(DY,SOUTH_NORTH_PATCH_END_STAG,coord_system=None):
    from iris import coords
    from numpy import arange
    Y_V=DY*(arange(0,SOUTH_NORTH_PATCH_END_STAG)-1)
    y_stag_coord=coords.AuxCoord(Y_V, standard_name='projection_y_coordinate', long_name='y', var_name='y', units='m', bounds=None, attributes=None, coord_system=coord_system)
    #y_stag_coord.add_dim_coord(south_north_stag,0)
    return y_stag_coord
    
def make_z_coordinate(filenames,**kwargs):
    from iris import coords
    z=calculate_wrf_geopotential_height(filenames,**kwargs)    
    z_coord=coords.AuxCoord(z.data, standard_name='geopotential_height', long_name='geopotential_height', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord

def make_z_xstag_coordinate(filenames,**kwargs):
    from iris import coords
    z=calculate_wrf_geopotential_height_xstag(filenames,**kwargs)
    z_coord=coords.AuxCoord(z, standard_name='geopotential_height', long_name='geopotential_height', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord
    
def make_z_ystag_coordinate(filenames,**kwargs):  
    from iris import coords
    z=calculate_wrf_geopotential_height_ystag(filenames,**kwargs)
    z_coord=coords.AuxCoord(z, standard_name='geopotential_height', long_name='geopotential_height', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord
    
    
def make_z_stag_coordinate(filenames,**kwargs):
    from iris import coords
    z=calculate_wrf_geopotential_height_stag(filenames,**kwargs)
    z_coord=coords.AuxCoord(z.data, standard_name='geopotential_height', long_name='z', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord
    
def make_p_coordinate(filenames,**kwargs):
    from iris import coords
    p=calculate_wrf_pressure(filenames,**kwargs)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord

def make_p_xstag_coordinate(filenames,**kwargs):
    from iris import coords
    p=calculate_wrf_pressure_xstag(filenames,**kwargs)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord
    
def make_p_ystag_coordinate(filenames,**kwargs):  
    from iris import coords
    p=calculate_wrf_pressure_ystag(filenames,**kwargs)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord
    
    
def make_p_stag_coordinate(filenames,**kwargs):
    from iris import coords
    p=calculate_wrf_pressure_stag(filenames,**kwargs)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord
    
def make_lon_coordinate(filenames,**kwargs):
    from iris import coords
    lon= loadwrfcube(filenames, 'XLONG',**kwargs)
    lon_coord=coords.AuxCoord(lon.data, standard_name=None, long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lon_coord
    
def make_lat_coordinate(filenames,**kwargs):
    from iris import coords
    lat= loadwrfcube(filenames, 'XLAT',**kwargs)
    lat_coord=coords.AuxCoord(lat.data, standard_name=None, long_name='latitude', var_name='latitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lat_coord
    
def make_lon_xstag_coordinate(filenames,**kwargs):
    from iris import coords
    lon= loadwrfcube(filenames, 'XLONG_U',**kwargs)
    lon_coord=coords.AuxCoord(lon.data, standard_name=None, long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lon_coord
    
def make_lat_xstag_coordinate(filenames,**kwargs):
    from iris import coords
    lat= loadwrfcube(filenames, 'XLAT_U',**kwargs)
    lat_coord=coords.AuxCoord(lat.data, standard_name=None, long_name='latitude', var_name='latitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lat_coord

def make_lon_ystag_coordinate(filenames,**kwargs):
    from iris import coords
    lon= loadwrfcube(filenames, 'XLONG_V',**kwargs)
    lon_coord=coords.AuxCoord(lon.data, standard_name=None, long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lon_coord
    
def make_lat_ystag_coordinate(filenames,**kwargs):
    from iris import coords
    lat= loadwrfcube(filenames, 'XLAT_V',**kwargs)
    lat_coord=coords.AuxCoord(lat.data, standard_name=None, long_name='latitude', var_name='latidude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return  lat_coord

