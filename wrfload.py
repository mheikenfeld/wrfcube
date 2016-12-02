from datetime import datetime

def loadwrfcube(filenames,variable,add_coordinates=None,slice_time=slice(None)):
#    print('loading ',variable)
    if type(filenames) is list:
        variable_cube=loadwrfcube_mult(filenames,variable,add_coordinates=add_coordinates,slice_time=slice_time)
    elif type(filenames) is str:
        variable_cube=loadwrfcube_single(filenames,variable,add_coordinates=add_coordinates,slice_time=slice_time)
    else:
        print('Type of input unknown: Must be str of list')
    variable_cube_data=variable_cube.data
    return variable_cube
    
def loadwrfcube_single(filenames,variable,add_coordinates=None,slice_time=slice(None)):
    from iris import load_cube    
    variable_cube=load_cube(filenames,variable)[slice_time]
    variable_cube=addcoordinates(filenames, variable,variable_cube,add_coordinates=addcoordinates, slice_time=slice_time)        
    return variable_cube
    
def loadwrfcube_mult(filenames,variable,add_coordinates=None,slice_time=slice(None)):
    from iris.cube import CubeList
    cube_list=[]
    for i in range(len(filenames)):
        cube_list.append(loadwrfcube_single(filenames[i],variable,add_coordinates=add_coordinates,slice_time=slice_time) )
    for member in cube_list:
        member.attributes={}
    variable_cubes=CubeList(cube_list)
    variable_cube=variable_cubes.concatenate_cube()
    return variable_cube

    
def loadwrfcube_dimcoord(filenames,variable,slice_time=slice(None)):
    from iris import load_cube     
    variable_cube=load_cube(filenames,variable)[slice_time]
    for coord in variable_cube.coords():
        variable_cube.remove_coord(coord.name())
    variable_cube=add_dim_coordinates(filenames, variable,variable_cube,slice_time=slice_time)
    variable_cube_data=variable_cube.data
    return variable_cube

    
def loadwrfcube_nocoord(filenames,variable,slice_time=slice(None)):
    from iris import load_cube   
    variable_cube=load_cube(filenames,variable)[slice_time]
    for coord in variable_cube.coords():
        variable_cube.remove_coord(coord.name())
    variable_cube_data=variable_cube.data
    return variable_cube
    

def derivewrfcube(filenames,variable,add_coordinates=None,slice_time=slice(None)):
    if type(filenames) is list:
        variable_cube=derivewrfcube_mult(filenames,variable,add_coordinates=add_coordinates,slice_time=slice_time)
    elif type(filenames) is str:
        variable_cube=derivewrfcube_single(filenames,variable,add_coordinates=add_coordinates,slice_time=slice_time)
    else:
        print('Type of input unknown: Must be str of list')
    return variable_cube

def derivewrfcube_mult(filenames,variable,add_coordinates=None,slice_time=slice(None)):
    from iris.cube import CubeList
    cube_list=[]
    for i in range(len(filenames)):
        cube_list.append(derivewrfcube_single(filenames[i],variable,add_coordinates=add_coordinates,slice_time=slice_time) )
    for member in cube_list:
        member.attributes={}
    variable_cubes=CubeList(cube_list)
    variable_cube=variable_cubes.concatenate_cube()
    variable_cube_data=variable_cube.data
    return variable_cube


def derivewrfcube_single(filenames,variable,add_coordinates=None,slice_time=slice(None)):
    if variable == 'potential temperature':
        variable_cube=calculate_wrf_potential_temperature(filenames,slice_time=slice_time)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'temperature':
        variable_cube=calculate_wrf_temperature(filenames,slice_time=slice_time)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'density':
        variable_cube=calculate_wrf_density(filenames,slice_time=slice_time)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'LWC':    
        variable_cube=calculate_wrf_LWC(filenames)
        #variable_cube=addcoordinates(filenames, 'QCLOUD',variable_cube,add_coordinates)
    elif variable == 'IWC':    
        variable_cube=calculate_wrf_IWC(filenames)
        #variable_cube=addcoordinates(filenames, 'QICE',variable_cube,add_coordinates)    
    elif variable == 'LWP':    
        variable_cube=calculate_wrf_LWP(filenames)
        #variable_cube=addcoordinates(filenames, 'OLR',variable_cube,add_coordinates)
    elif variable == 'IWP':    
        variable_cube=calculate_wrf_IWP(filenames)
        #variable_cube=addcoordinates(filenames, 'OLR',variable_cube,add_coordinates)
    elif variable == 'IWV':    
        variable_cube=calculate_wrf_IWV(filenames)
        #variable_cube=addcoordinates(filenames, 'OLR',variable_cube,add_coordinates)
    elif variable == 'airmass':    
        variable_cube=calculate_wrf_airmass(filenames)
    elif variable == 'geopotential':    
        variable_cube=calculate_wrf_geopotential(filenames,slice_time=slice_time)
        variable_cube=remove_all_coordinates(variable_cube)
        variable_cube=addcoordinates(filenames, 'QICE',variable_cube,add_coordinates=add_coordinates,slice_time=slice_time)
    elif variable == 'pressure':    
        variable_cube=calculate_wrf_pressure(filenames,slice_time=slice_time)
        #variable_cube=addcoordinates(filenames, 'T',variable_cube,add_coordinates=add_coordinates,slice_time=slice_time)
        #print(variable_cube)
    elif variable == 'relative humidity':    
        variable_cube=calculate_wrf_relativehumidity(filenames)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'w_at_T':    
        variable_cube=calculate_wrf_w_at_T(filenames)
        variable_cube=addcoordinates(filenames, 'T',variable_cube,add_coordinates=add_coordinates,slice_time=slice_time)
    elif variable == 'surface precipitation':
        variable_cube=calculate_wrf_surface_precipitation(filenames)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'maximum reflectivity':    
        variable_cube=calculate_wrf_maximum_reflectivity(filenames)
    else:
        raise NameError(variable, 'is not a known variable') 
    return variable_cube
    
def calculate_wrf_surface_precipitation(filenames):
    import numpy as np
    RAINNC= loadwrfcube(filenames, 'RAINNC')
    dt=(RAINNC.coords('Time')[0].points[1]-RAINNC.coords('Time')[0].points[0])*24
    rainnc_inst=np.concatenate((RAINNC.data[[1],:,:]-RAINNC.data[[0],:,:],RAINNC.data[1:,:,:]-RAINNC.data[0:-1:,:,:]),axis=0)/dt
    RAINNC_inst=RAINNC
    RAINNC_inst.data=rainnc_inst
    RAINNC_inst.rename('surface precipitation')
    RAINNC_inst.units= 'mm/h'
    return RAINNC_inst

def variable_list(filenames):
    from netCDF4 import Dataset
    variable_list = list(Dataset(filenames).variables)
    return variable_list

def calculate_wrf_potential_temperature(filenames,slice_time=slice(None)):
    from iris import coords
    T= loadwrfcube(filenames, 'T',slice_time=slice_time)
    T0 = coords.AuxCoord(300.0,long_name='reference_Temperature', units='K')
    theta=T+T0;
    theta.rename('potential temperature')
    return theta
    
def calculate_wrf_temperature(filenames,slice_time=slice(None)):
    from iris import coords, cube
    theta= calculate_wrf_potential_temperature(filenames,slice_time=slice_time)  
    p = derivewrfcube(filenames,'pressure',slice_time=slice_time)
    p0 =coords.AuxCoord(1000.0,long_name='reference_pressure', units='hPa')
    p0.convert_units(p.units)
    p1=p/p0
    exp=(287.05 / 1005.0)
    T=theta*(p1**exp) #work around iris issue here by loading one of the cubes into numpy array..
    T.rename('air_temperature')
    return T
    
def calculate_wrf_relativehumidity(filenames,slice_time=slice(None)):
    from iris import cube
    QVAPOR=loadwrfcube(filenames, 'QVAPOR',slice_time=slice_time).data
    T=calculate_wrf_temperature(filenames,slice_time=slice_time).data
    p=calculate_wrf_pressure(filenames,slice_time=slice_time)
    p.convert_units('Pa')
    p=p.data
    rh=calculate_RH(QVAPOR,T,p)
    RH=cube.Cube(rh, units='percent',long_name='realtive humidity')
    return RH   
    
def calculate_RH(QVAPOR,T,p):
    from numpy import exp,maximum,minimum
    ES=1e2*6.1094*exp(17.625*(T-273.15)/(T-273.15+243.04))
    QVS = 0.622*ES/ (p- (1.0-0.622)*ES)
    RH = 100*maximum(minimum(QVAPOR/QVS,1.0),0.0)
    return RH   
    
def calculate_wrf_LWC(filenames):
    QCLOUD=loadwrfcube(filenames, 'QCLOUD')
    QRAIN=loadwrfcube(filenames, 'QRAIN')
    LWC=QCLOUD+QRAIN
    LWC.rename('liquid water content')
    return LWC   
#    
def calculate_wrf_IWC(filenames):    
    QICE=loadwrfcube(filenames, 'QICE')
    QSNOW=loadwrfcube(filenames, 'QSNOW')
    QGRAUP=loadwrfcube(filenames, 'QGRAUP')
    IWC=QICE+QSNOW+QGRAUP
    IWC.rename('Ice water content')
    return IWC
    
def calculate_wrf_airmass(filenames):
    rho=calculate_wrf_density(filenames)
    z_h=calculate_wrf_geopotential_stag(filenames)
    Airmass=rho*(z_h[:,1:,:,:].data-z_h[:,0:-1,:,:].data)
    Airmass.rename('mass of air')
    return Airmass
    
def calculate_wrf_LWP(filenames):
    from iris.analysis import SUM
    LWC=calculate_wrf_LWC(filenames)
    Airmass=calculate_wrf_airmass(filenames)
    LWP=(LWC*Airmass).collapsed(('bottom_top'),SUM)
    LWP.rename('liquid water path')
    return LWP   
#    
def calculate_wrf_IWP(filenames):    
    from iris.analysis import SUM
    IWC=calculate_wrf_IWC(filenames)
    Airmass=calculate_wrf_airmass(filenames)
    IWP=(IWC*Airmass).collapsed(('bottom_top'),SUM)
    IWP.rename('ice water path')
    return IWP
    
def calculate_wrf_IWV(filenames):    
    from iris.analysis import SUM
    QVAPOR=loadwrfcube(filenames,'QVAPOR')
    Airmass=calculate_wrf_airmass(filenames)
    IWV=(QVAPOR*Airmass).collapsed(('bottom_top'),SUM)
    IWV.rename('integrated water vapor')
    return IWV
    
def calculate_wrf_maximum_reflectivity(filenames):
    from iris.analysis import MAX
    REFL_10CM=loadwrfcube(filenames,'REFL_10CM')
    MAX_REFL_10CM=REFL_10CM.collapsed('bottom_top', MAX)
    MAX_REFL_10CM.rename('maximum reflectivity')
    return MAX_REFL_10CM
    
def calculate_wrf_w_at_T(filenames):
    from iris import cube
    w=loadwrfcube(filenames, 'W')
    w_at_T = cube.Cube(0.5*(w[:,:-1,:,:].data+w[:,1:,:,:].data),var_name='w',long_name='vertical velocity on T grid', units='m/s')
    return w_at_T

def calculate_wrf_density(filenames,slice_time=slice(None)):
    from iris import coords
    if ('ALT' in variable_list(filenames)):
        alt=loadwrfcube(filenames,'ALT',slice_time=slice_time)
        rho=alt**(-1)
    else: 
       R=coords.AuxCoord(287.058,long_name='Specific gas constant for air',units='Joule kg^-1 K^-1')
       T=calculate_wrf_temperature(filenames,slice_time=slice_time)
       p=derivewrfcube(filenames,'pressure',slice_time=slice_time)
       rho=p*((R*T)**-1)
       rho.rename('air_density')
    return rho
#    
def calculate_wrf_pressure(filenames,slice_time=slice(None)):
    P= loadwrfcube(filenames, 'P',slice_time=slice_time)
    PB= loadwrfcube(filenames, 'PB',slice_time=slice_time)
    p=P + PB 
    p.rename('pressure')
    return p
    
    #    
def calculate_wrf_pressure_xstag(filenames):
    p=calculate_wrf_pressure(filenames)
    p_xstag = 0.5*(p[:,:,:-1,:]+p[:,:,1:,:])
    p_xstag.rename('pressure')
    return p

#    
def calculate_wrf_pressure_ystag(filenames):
    p=calculate_wrf_pressure(filenames)
    p_ystag = 0.5*(p[:,:,:-1,:]+p[:,:,1:,:])
    p_ystag.rename('pressure')
    return p_ystag

    
def calculate_wrf_pressure_stag(filenames,slice_time=slice(None)):
    PH= loadwrfcube(filenames, 'PH',slice_time=slice_time)
    PHB= loadwrfcube(filenames,'PHB',slice_time=slice_time)
    pH=PH + PHB
    pH.rename('pressure')
    return pH
    
    
    
def calculate_wrf_geopotential_stag(filenames,slice_time=slice(None)):
    from iris import coords
    pH=calculate_wrf_pressure_stag(filenames,slice_time=slice_time)
    g = coords.AuxCoord(9.81,long_name='acceleration', units='m s^-2')
    zH=pH/g
    zH.rename('geopotential')
    return zH

def calculate_wrf_geopotential(filenames,slice_time=slice(None)):
    zH=calculate_wrf_geopotential_stag(filenames,slice_time=slice_time)
    z = 0.5*(zH[:,:-1,:,:]+zH.data[:,1:,:,:])
    z.rename('geopotential')
    return z
    
def calculate_wrf_geopotential_ystag(filenames):
    z=calculate_wrf_geopotential(filenames)
    z_ystag = 0.5*(z[:,:,:-1,:]+z[:,:,1:,:])
    z_ystag.rename('geopotential')
    return z_ystag
    
def calculate_wrf_geopotential_xstag(filenames):
    z=calculate_wrf_geopotential(filenames)
    z_xstag = 0.5*(z[:,:,:,:-1]+z[:,:,:,1:])
    z_xstag.rename('geopotential')
    return z_xstag

def remove_all_coordinates(variable_cube):
    for coordinate in variable_cube.coords():
        variable_cube.remove_coord(coordinate.name())
    return variable_cube    
    
    
def addcoordinates(filenames, variable,variable_cube,add_coordinates=None,slice_time=slice(None)):
    if add_coordinates==None:
        variable_cube=add_dim_coordinates(filenames, variable,variable_cube,slice_time=slice_time)
    else:
        variable_cube=add_dim_coordinates(filenames, variable,variable_cube,slice_time=slice_time)
        variable_cube=add_aux_coordinates(filenames, variable,variable_cube,add_coordinates)
    return variable_cube

def add_dim_coordinates(filenames, variable,variable_cube,slice_time=slice(None)):
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
           time=make_time_coord(filenames,slice_time=slice_time)
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
           west_east_stag=make_westeast_coord(DX,WEST_EAST_PATCH_END_STAG)
           variable_cube.add_dim_coord(west_east_stag,dim)
        elif variable_dimensions[dim]=='south_north_stag':
           south_north_stag=make_southnorth_coord(DY, SOUTH_NORTH_PATCH_END_STAG)
           variable_cube.add_dim_coord(south_north_stag,dim)
        elif variable_dimensions[dim]=='bottom_top_stag':
           bottom_top_stag=make_bottom_top_stag_coordinate(BOTTOM_TOP_PATCH_END_STAG)   
           variable_cube.add_dim_coord(bottom_top_stag,dim)
    return variable_cube        
    
def add_aux_coordinates(filenames, variable,variable_cube,add_coordinates):
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
    coords=variable_cube.coords()
    if type(add_coordinates)!=list:
        add_coordinates1=add_coordinates
        add_coordinates=[]
        add_coordinates.append(add_coordinates1)
    for coordinate in add_coordinates:
        if coordinate=='xy':
            for dim in range(len(coords)):
                if (coords[dim].name()=='west_east'):
                    x_coord=make_x_coord(DX,WEST_EAST_PATCH_END_UNSTAG)
                    variable_cube.add_aux_coord(x_coord,dim)
                elif (coords[dim].name()=='south_north'):
                    y_coord=make_y_coord(DY, SOUTH_NORTH_PATCH_END_UNSTAG)
                    variable_cube.add_aux_coord(y_coord,dim)
                elif (coords[dim].name()=='west_east_stag'):
                    x_stag_coord=make_x_stag_coord(DX,WEST_EAST_PATCH_END_STAG)
                    variable_cube.add_aux_coord(x_stag_coord,dim)
                elif coords[dim].name()=='south_north_stag':
                    y_stag_coord=make_y_stag_coord(DY, SOUTH_NORTH_PATCH_END_STAG)
                    variable_cube.add_aux_coord(y_stag_coord,dim)
                    
        if coordinate=='z':    
            if (coords[0].name()=='Time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                z_coord=make_z_coordinate(filenames,coords)
                variable_cube.add_aux_coord(z_coord,(0,1,2,3))
            elif (coords[0].name()=='Time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east_stag'):
                z_coord=make_z_xstag_coordinate(filenames)     
                variable_cube.add_aux_coord(z_coord)
            elif (coords[0].name()=='Time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north_stag' and coords[3].name()=='west_east'):
                z_coord=make_z_ystag_coordinate(filenames)
                variable_cube.add_aux_coord(z_coord,(0,1,2,3))
            elif (coords[0].name()=='Time' and coords[1].name()=='bottom_top_stag' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                z_stag_coord=make_z_stag_coordinate(filenames)   
                variable_cube.add_aux_coord(z_stag_coord,(0,1,2,3))
            else:
                print("no z coordinates added")
                
        if coordinate=='pressure' :      
            if (coords[0].name()=='Time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                p_coord=make_p_coordinate(filenames,coords)
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            elif (coords[0].name()=='Time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east_stag'):
                p_coord=make_p_xstag_coordinate(filenames)
                variable_cube.add_aux_coord(p_coord)
            elif (coords[0].name()=='Time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north_stag' and coords[3].name()=='west_east'):
                p_coord=make_p_ystag_coordinate(filenames)
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            elif (coords[0].name()=='Time' and coords[1].name()=='bottom_top_stag' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                p_coord=make_p_stag_coordinate(filenames)   
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            else:
                print("p coordinates added")
                
        if (coordinate=='zp' or coordinate=='pz'):    
            if (coords[0].name()=='Time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                z_coord=make_z_coordinate(filenames,coords)
                variable_cube.add_aux_coord(z_coord,(0,1,2,3))
                p_coord=make_p_coordinate(filenames,coords)
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            elif (coords[0].name()=='Time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north' and coords[3].name()=='west_east_stag'):
                z_coord=make_z_xstag_coordinate(filenames)     
                variable_cube.add_aux_coord(z_coord)
                p_coord=make_p_xstag_coordinate(filenames)
                variable_cube.add_aux_coord(p_coord)
            elif (coords[0].name()=='Time' and coords[1].name()=='bottom_top' and coords[2].name()=='south_north_stag' and coords[3].name()=='west_east'):
                z_coord=make_z_ystag_coordinate(filenames)
                variable_cube.add_aux_coord(z_coord,(0,1,2,3))
                p_coord=make_p_ystag_coordinate(filenames)
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            elif (coords[0].name()=='Time' and coords[1].name()=='bottom_top_stag' and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                z_stag_coord=make_z_stag_coordinate(filenames)   
                variable_cube.add_aux_coord(z_stag_coord,(0,1,2,3))
                p_coord=make_p_stag_coordinate(filenames)   
                variable_cube.add_aux_coord(p_coord,(0,1,2,3))
            else:
                print("no z and p coordinates added")
                
        if coordinate=='latlon':    
            if (coords[0].name()=='Time' and (coords[1].name()=='bottom_top' or 'bottom_top_stag') and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                lat_coord=make_lat_coordinate(filenames)
                lon_coord=make_lon_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,2,3))
                variable_cube.add_aux_coord(lon_coord,(0,2,3))            
            elif (coords[0].name()=='Time' and (coords[1].name()=='bottom_top' or 'bottom_top_stag') and coords[2].name()=='south_north' and coords[3].name()=='west_east_stag'):
                lat_coord=make_lat_xstag_coordinate(filenames)
                lon_coord=make_lon_xstag_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,2,3))
                variable_cube.add_aux_coord(lon_coord,(0,2,3))
            elif (coords[0].name()=='Time' and (coords[1].name()=='bottom_top' or 'bottom_top_stag') and coords[2].name()=='south_north_stag' and coords[3].name()=='west_east'):
                lat_coord=make_lat_ystag_coordinate(filenames)
                lon_coord=make_lon_ystag_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,2,3))
                variable_cube.add_aux_coord(lon_coord,(0,2,3))
            elif (coords[0].name()=='Time'  and coords[1].name()=='south_north' and coords[2].name()=='west_east'):
                lat_coord=make_lat_coordinate(filenames)
                lon_coord=make_lon_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,1,2))
                variable_cube.add_aux_coord(lon_coord,(0,1,2))            
            elif (coords[0].name()=='Time'  and coords[1].name()=='south_north' and coords[2].name()=='west_east_stag'):
                lat_coord=make_lat_xstag_coordinate(filenames)
                lon_coord=make_lon_xstag_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,1,2))
                variable_cube.add_aux_coord(lon_coord,(0,1,2))
            elif (coords[0].name()=='Time' and coords[1].name()=='south_north_stag' and coords[2].name()=='west_east'):
                lat_coord=make_lat_ystag_coordinate(filenames)
                lon_coord=make_lon_ystag_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,1,2))
                variable_cube.add_aux_coord(lon_coord,(0,1,2))
            else:
                print("no lat/lon coordinates added")
    return variable_cube
    
def make_time_coord(filenames,slice_time=slice(None)):
    from iris import load_cube,coords
    from datetime import datetime,timedelta
    from numpy import empty
    Times= load_cube(filenames, 'Times')
    filetimes = Times.data[slice_time]
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
    for i in range(len(timeobjlist)):
        time_days[i]=(timeobjlist[i] - datetime(1970,1,1)).total_seconds() / timedelta(1).total_seconds()
    time_coord=coords.DimCoord(time_days, standard_name=None, long_name='Time', var_name='Time', units='days since 1970-01-01', bounds=None, attributes=None, coord_system=None, circular=False)
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
    west_east_stag=coords.DimCoord(WEST_EAST_U, standard_name=None, long_name='west_east', var_name='west_east', units='1', bounds=None, attributes=None, coord_system=None, circular=False)
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
    south_north_stag=coords.DimCoord(SOUTH_NORTH_V, standard_name=None, long_name='south_north', var_name='south_north', units='1', bounds=None, attributes=None, coord_system=None, circular=False)
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

def make_x_coord(DX,WEST_EAST_PATCH_END_UNSTAG):
    from iris import coords
    from numpy import arange
    X=DX*(arange(0,WEST_EAST_PATCH_END_UNSTAG)-0.5)
    x_coord=coords.AuxCoord(X, standard_name=None, long_name='x', var_name='x', units='km', bounds=None, attributes=None, coord_system=None)
    #x_coord.add_dim_coord(west_east,0)
    return x_coord

def make_x_stag_coord(DX,WEST_EAST_PATCH_END_STAG):
    from iris import coords
    from numpy import arange
    X_U=DX*(arange(0,WEST_EAST_PATCH_END_STAG)-1)
    x_stag_coord=coords.AuxCoord(X_U, standard_name=None, long_name='x', var_name='x', units='km', bounds=None, attributes=None, coord_system=None)
    #x_stag_coord.add_dim_coord(west_east_stag,0)
    return x_stag_coord

def make_y_coord(DY,SOUTH_NORTH_PATCH_END_UNSTAG):
    from iris import coords
    from numpy import arange    #DY=attributes['DY']
    Y=DY*(arange(0,SOUTH_NORTH_PATCH_END_UNSTAG)-0.5)
    y_coord=coords.AuxCoord(Y, standard_name=None, long_name='y', var_name='y', units='km', bounds=None, attributes=None, coord_system=None)
    #y_coord.add_dim_coord(south_north,0)
    return y_coord
    
def make_y_stag_coord(DY,SOUTH_NORTH_PATCH_END_STAG):
    from iris import coords
    from numpy import arange
    Y_V=DY*(arange(0,SOUTH_NORTH_PATCH_END_STAG)-1)
    y_stag_coord=coords.AuxCoord(Y_V, standard_name=None, long_name='y', var_name='y', units='km', bounds=None, attributes=None, coord_system=None)
    #y_stag_coord.add_dim_coord(south_north_stag,0)
    return y_stag_coord
    
def make_z_coordinate(filenames,coordinates):
    from iris import coords
    z=calculate_wrf_geopotential(filenames)    
    z_coord=coords.AuxCoord(z.data, standard_name=None, long_name='geopotential', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord

def make_z_xstag_coordinate(filenames):
    from iris import coords
    z=calculate_wrf_geopotential_xstag(filenames)
    z_coord=coords.AuxCoord(z.data, standard_name=None, long_name='z', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord
    
def make_z_ystag_coordinate(filenames):  
    from iris import coords
    z=calculate_wrf_geopotential_ystag(filenames)
    z_coord=coords.AuxCoord(z.data, standard_name=None, long_name='z', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord
    
    
def make_z_stag_coordinate(filenames):
    from iris import coords
    z=calculate_wrf_geopotential_stag(filenames)
    z_coord=coords.AuxCoord(z.data, standard_name=None, long_name='z', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord
    
def make_p_coordinate(filenames,coordinates):
    from iris import coords
    p=calculate_wrf_pressure(filenames)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord

def make_p_xstag_coordinate(filenames):
    from iris import coords
    p=calculate_wrf_pressure_xstag(filenames)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord
    
def make_p_ystag_coordinate(filenames):  
    from iris import coords
    p=calculate_wrf_pressure_ystag(filenames)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord
    
    
def make_p_stag_coordinate(filenames):
    from iris import coords
    p=calculate_wrf_pressure_stag(filenames)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord
    
def make_lon_coordinate(filenames):
    from iris import coords
    lon= loadwrfcube(filenames, 'XLONG')
    lon_coord=coords.AuxCoord(lon.data, standard_name=None, long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lon_coord
    
def make_lat_coordinate(filenames):
    from iris import coords
    lat= loadwrfcube(filenames, 'XLAT')
    lat_coord=coords.AuxCoord(lat.data, standard_name=None, long_name='latitude', var_name='latitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lat_coord
    
def make_lon_xstag_coordinate(filenames):
    from iris import coords
    lon= loadwrfcube(filenames, 'XLONG_U')
    lon_coord=coords.AuxCoord(lon.data, standard_name=None, long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lon_coord
    
def make_lat_xstag_coordinate(filenames):
    from iris import coords
    lat= loadwrfcube(filenames, 'XLAT_U')
    lat_coord=coords.AuxCoord(lat.data, standard_name=None, long_name='latitude', var_name='latitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lat_coord

def make_lon_ystag_coordinate(filenames):
    from iris import coords
    lon= loadwrfcube(filenames, 'XLONG_V')
    lon_coord=coords.AuxCoord(lon.data, standard_name=None, long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lon_coord
    
def make_lat_ystag_coordinate(filenames):
    from iris import coords
    lat= loadwrfcube(filenames, 'XLAT_V')
    lat_coord=coords.AuxCoord(lat.data, standard_name=None, long_name='latitude', var_name='latidude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return  lat_coord

