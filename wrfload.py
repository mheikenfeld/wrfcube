def loadwrfcube(filenames,variable,add_coordinates=None):
    if type(filenames) is list:
        variable_cube=loadwrfcube_mult(filenames,variable,add_coordinates)
    elif type(filenames) is str:
        variable_cube=loadwrfcube_single(filenames,variable,add_coordinates)
    else:
        print('Type of input unknown: Must be str of list')
    return variable_cube
    
def loadwrfcube_single(filenames,variable,add_coordinates=None):
    from iris import load 
    variable_cube=load(filenames,variable)[0]
    variable_cube=addcoordinates(filenames, variable,variable_cube,add_coordinates)        
    return variable_cube
    
def loadwrfcube_mult(filenames,variable,add_coordinates=None):
    from iris.cube import CubeList
    cube_list=[]
    for i in range(len(filenames)):
        cube_list.append(loadwrfcube_single(filenames[i],variable,add_coordinates) )
    for member in cube_list:
        member.attributes={}
    variable_cubes=CubeList(cube_list)
    variable_cube=variable_cubes.concatenate_cube()
    return variable_cube

    
def loadwrfcube_dimcoord(filenames,variable):
    from iris import load     
    variable_cube=load(filenames,variable)[0]
    for coord in variable_cube.coords():
        variable_cube.remove_coord(coord.name())
    variable_cube=add_dim_coordinates(filenames, variable,variable_cube)
    return variable_cube

    
def loadwrfcube_nocoord(filenames,variable):
    from iris import load     
    variable_cube=load(filenames,variable)[0]
    for coord in variable_cube.coords():
        variable_cube.remove_coord(coord.name())
    return variable_cube
    

def derivewrfcube(filenames,variable,add_coordinates=None):
    if type(filenames) is list:
        variable_cube=derivewrfcube_mult(filenames,variable,add_coordinates)
    elif type(filenames) is str:
        variable_cube=derivewrfcube_single(filenames,variable,add_coordinates)
    else:
        print('Type of input unknown: Must be str of list')
    return variable_cube

def derivewrfcube_mult(filenames,variable,add_coordinates=None):
    from iris.cube import CubeList
    cube_list=[]
    for i in range(len(filenames)):
        cube_list.append(derivewrfcube_single(filenames[i],variable,add_coordinates) )
    for member in cube_list:
        member.attributes={}
    variable_cubes=CubeList(cube_list)
    variable_cube=variable_cubes.concatenate_cube()
    return variable_cube


def derivewrfcube_single(filenames,variable,add_coordinates=None):
    from iris import load,cube
    if variable == 'potential temperature':
        variable_cube=calculate_wrf_potential_temperature(filenames)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'temperature':
        variable_cube=calculate_wrf_temperature(filenames)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'density':
        variable_cube=calculate_wrf_density(filenames)
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
    elif variable == 'geopotential':    
        variable_cube=calculate_wrf_geopotential(filenames)
        variable_cube=addcoordinates(filenames, 'QICE',variable_cube,add_coordinates)
    elif variable == 'pressure':    
        variable_cube=calculate_wrf_pressure(filenames)
        variable_cube=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'relative humidity':    
        variable_cube=calculate_wrf_relativehumidity(filenames)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
    elif variable == 'w_at_T':    
        variable_cube=calculate_wrf_w_at_T(filenames)
        #variable_cube_out=addcoordinates(filenames, 'T',variable_cube,add_coordinates)
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


def calculate_wrf_potential_temperature(filename):
    from iris import coords
    THETA= loadwrfcube(filename, 'T')
    T0 = coords.AuxCoord(300.0,long_name='reference_Temperature',units='K')
    theta=THETA+T0;
    theta.rename('Potential temperature')
    theta.attributes=THETA.attributes
    return theta
    
def calculate_wrf_temperature(filename):
    from iris import coords
    theta= calculate_wrf_potential_temperature(filename)
    p = calculate_wrf_pressure(filename)
    p0 =coords.AuxCoord(1000.0,long_name='reference_pressure', units='hPa')
    p0.convert_units(p.units)
    p1=p/p0
    T=theta
    T.data=theta.data*(p1.data**(287.05 / 1005))
    T.rename('temperature')
    T.attributes=theta.attributes
    return T
    
def calculate_wrf_relativehumidity(filename):
    from iris import cube
    QVAPOR=loadwrfcube(filename, 'QVAPOR').data
    T=calculate_wrf_temperature(filename).data
    p=calculate_wrf_pressure(filename)
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
    
def calculate_wrf_LWC(filename):
    from iris import load
    QCLOUD=loadwrfcube(filename, 'QCLOUD')
    QRAIN=loadwrfcube(filename, 'QRAIN')
    LWC=QCLOUD+QRAIN
    LWC.rename('Liquid water content')
    return LWC   
#    
def calculate_wrf_IWC(filename):    
    QICE=loadwrfcube(filename, 'QICE')
    QSNOW=loadwrfcube(filename, 'QSNOW')
    QGRAUP=loadwrfcube(filename, 'QGRAUP')
    IWC=QICE+QSNOW+QGRAUP
    IWC.rename('Ice water content')
    return IWC
    
def calculate_wrf_LWP(filename):
    from iris import cube, coords    
    import numpy as np    
    LWC=calculate_wrf_LWC(filename)
    rho=calculate_wrf_density(filename)
    z_h=calculate_wrf_geopotential_stag(filename)
    LWP=cube.Cube(np.sum(((LWC)*rho*(z_h[:,1:,:,:]-z_h[:,0:-1,:,:])).data,axis=1),long_name='liquid water path', var_name='liquid_water_path', units='kg m^-2')
    lat_coord=coords.AuxCoord(LWC.coords('latitude')[0].points, standard_name='latitude', long_name='latitude', var_name='latitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    lon_coord=coords.AuxCoord(LWC.coords('longitude')[0].points, standard_name='longitude', long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None) 
    LWP.add_aux_coord(lat_coord,(0,1,2))
    LWP.add_aux_coord(lon_coord,(0,1,2))            

    return LWP   
#    
def calculate_wrf_IWP(filename):    
    from iris import cube, coords    
    import numpy as np    
    IWC=calculate_wrf_IWC(filename)
    rho=calculate_wrf_density(filename)
    z_h=calculate_wrf_geopotential_stag(filename)
    IWP=cube.Cube(np.sum((IWC*rho*(z_h[:,1:,:,:]-z_h[:,0:-1,:,:])).data,axis=1),long_name='ice water path', var_name='ice_water_path', units='kg m^-2')
    lat_coord=coords.AuxCoord(IWC.coords('latitude')[0].points, standard_name='latitude', long_name='latitude', var_name='latitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    lon_coord=coords.AuxCoord(IWC.coords('longitude')[0].points, standard_name='longitude', long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None) 
    IWP.add_aux_coord(lat_coord,(0,1,2))
    IWP.add_aux_coord(lon_coord,(0,1,2))            

    return IWP
    
def calculate_wrf_maximum_reflectivity(filename):
    from iris.analysis import MAX
    REFL_10CM=loadwrfcube(filename,'REFL_10CM')
    MAX_REFL_10CM=REFL_10CM.collapsed('bottom_top', MAX)
    MAX_REFL_10CM.rename('maximum reflectivity')
    return MAX_REFL_10CM
    
def calculate_wrf_w_at_T(filename):
    w=loadwrfcube(filename, 'W')
    w_at_T = 0.5*(w[:,:-1,:,:]+w[:,1:,:,:])
    return w_at_T

def calculate_wrf_density(filename):    
    T=calculate_wrf_temperature(filename)
    p=calculate_wrf_pressure(filename)
    rho=p/ (287.*T)
    return rho
#    
def calculate_wrf_pressure(filename):
    P= loadwrfcube_nocoord(filename, 'P')
    PB= loadwrfcube_nocoord(filename, 'PB')
    p=P + PB
    p.rename('pressure')
    return p
    
    #    
def calculate_wrf_pressure_xstag(filename):
    p=calculate_wrf_pressure(filename)
    p_xstag = 0.5*(p[:,:,:-1,:]+p[:,:,1:,:])
    p_xstag.rename('pressure')
    return p

#    
def calculate_wrf_pressure_ystag(filename):
    p=calculate_wrf_pressure(filename)
    p_ystag = 0.5*(p[:,:,:-1,:]+p[:,:,1:,:])
    p_ystag.rename('pressure')
    return p_ystag

    
def calculate_wrf_pressure_stag(filename):
    PH= loadwrfcube_nocoord(filename, 'PH')
    PHB= loadwrfcube_nocoord(filename,'PHB')
    pH=PH + PHB
    pH.rename('pressure')
    return pH
    
    
    
def calculate_wrf_geopotential_stag(filename):
    from iris import coords
    pH=calculate_wrf_pressure_stag(filename)
    g = coords.AuxCoord(9.81,long_name='acceleration', units='m s^-2')
    zH=pH/g
    zH.rename('geopotential')
    return zH

def calculate_wrf_geopotential(filename):
    zH=calculate_wrf_geopotential_stag(filename)
    z = 0.5*(zH[:,:-1,:,:]+zH[:,1:,:,:])
    z.rename('geopotential')
    return z
    
def calculate_wrf_geopotential_ystag(filename):
    z=calculate_wrf_geopotential(filename)
    z_ystag = 0.5*(z[:,:,:-1,:]+z[:,:,1:,:])
    z_ystag.rename('geopotential')
    return z_ystag
    
def calculate_wrf_geopotential_xstag(filename):
    z=calculate_wrf_geopotential(filename)
    z_xstag = 0.5*(z[:,:,:,:-1]+z[:,:,:,1:])
    z_xstag.rename('geopotential')
    return z_xstag

def addcoordinates(filenames, variable,variable_cube,add_coordinates=None):
    if add_coordinates==None:
        variable_cube=add_dim_coordinates(filenames, variable,variable_cube)
    else:
        variable_cube=add_dim_coordinates(filenames, variable,variable_cube)
        variable_cube=add_aux_coordinates(filenames, variable,variable_cube,add_coordinates)
    return variable_cube

def add_dim_coordinates(filenames, variable,variable_cube):
    from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
    from iris import load
    variable_cube_dim= load(filenames, variable)[0]
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
    from iris import load
    variable_cube_dim= load(filenames, variable)[0]
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
    for coordinate in add_coordinates:
        
        if add_coordinates=='xy':
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
                    
        if add_coordinates=='z':    
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
                
        if add_coordinates=='pressure' :      
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
                
        if add_coordinates=='latlon':    
            if (coords[0].name()=='Time' and (coords[1].name()=='bottom_top' or 'bottom_top_stag') and coords[2].name()=='south_north' and coords[3].name()=='west_east'):
                lat_coord=make_lat_coordinate(filenames)
                lon_coord=make_lon_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,2,3))
                variable_cube.add_aux_coord(lon_coord,(0,2,3))            
            elif (coords[0].name()=='Time' and (coords[1].name()=='bottom_top' or 'bottom_top_stag') and coords[2].name()=='south_north' and coords[3].name()=='west_east_stag'):
                lat_coord=make_lat_stagx_coordinate(filenames)
                lon_coord=make_lon_stagx_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,2,3))
                variable_cube.add_aux_coord(lon_coord,(0,2,3))
            elif (coords[0].name()=='Time' and (coords[1].name()=='bottom_top' or 'bottom_top_stag') and coords[2].name()=='south_north_stag' and coords[3].name()=='west_east'):
                lat_coord=make_lat_stagy_coordinate(filenames)
                lon_coord=make_lon_stagy_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,2,3))
                variable_cube.add_aux_coord(lon_coord,(0,2,3))
            elif (coords[0].name()=='Time'  and coords[1].name()=='south_north' and coords[2].name()=='west_east'):
                lat_coord=make_lat_coordinate(filenames)
                lon_coord=make_lon_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,1,2))
                variable_cube.add_aux_coord(lon_coord,(0,1,2))            
            elif (coords[0].name()=='Time'  and coords[1].name()=='south_north' and coords[2].name()=='west_east_stag'):
                lat_coord=make_lat_stagx_coordinate(filenames)
                lon_coord=make_lon_stagx_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,1,2))
                variable_cube.add_aux_coord(lon_coord,(0,1,2))
            elif (coords[0].name()=='Time' and coords[1].name()=='south_north_stag' and coords[2].name()=='west_east'):
                lat_coord=make_lat_stagy_coordinate(filenames)
                lon_coord=make_lon_stagy_coordinate(filenames)   
                variable_cube.add_aux_coord(lat_coord,(0,1,2))
                variable_cube.add_aux_coord(lon_coord,(0,1,2))
            else:
                print("no lat/lon coordinates added")
    return variable_cube
    
def make_time_coord(filename):
    from iris import load,coords
    from datetime import datetime,timedelta
    from numpy import empty
    Times= load(filename, 'Times')
    filetimes = Times[0].data
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
    
def make_z_coordinate(filename,coordinates):
    from iris import coords
    z=calculate_wrf_geopotential(filename)    
    z_coord=coords.AuxCoord(z.data, standard_name=None, long_name='geopotential', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord

def make_z_xstag_coordinate(filename):
    from iris import coords
    z=calculate_wrf_geopotential_xstag(filename)
    z_coord=coords.AuxCoord(z.data, standard_name=None, long_name='z', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord
    
def make_z_ystag_coordinate(filename):  
    from iris import coords
    z=calculate_wrf_geopotential_ystag(filename)
    z_coord=coords.AuxCoord(z.data, standard_name=None, long_name='z', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord
    
    
def make_z_stag_coordinate(filename):
    from iris import coords
    z=calculate_wrf_geopotential_stag(filename)
    z_coord=coords.AuxCoord(z.data, standard_name=None, long_name='z', var_name='z', units='m', bounds=None, attributes=None, coord_system=None)
    return z_coord
    
def make_p_coordinate(filename,coordinates):
    from iris import coords
    p=calculate_wrf_pressure(filename)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord

def make_p_xstag_coordinate(filename):
    from iris import coords
    p=calculate_wrf_pressure_xstag(filename)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord
    
def make_p_ystag_coordinate(filename):  
    from iris import coords
    p=calculate_wrf_pressure_ystag(filename)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord
    
    
def make_p_stag_coordinate(filename):
    from iris import coords
    p=calculate_wrf_pressure_stag(filename)
    p_coord=coords.AuxCoord(p.data, standard_name=None, long_name='pressure', var_name='pressure', units='Pa', bounds=None, attributes=None, coord_system=None)
    return p_coord
    
def make_lon_coordinate(filename):
    from iris import coords
    lon= loadwrfcube(filename, 'XLONG')
    lon_coord=coords.AuxCoord(lon.data, standard_name=None, long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lon_coord
    
def make_lat_coordinate(filename):
    from iris import coords
    lat= loadwrfcube(filename, 'XLAT')
    lat_coord=coords.AuxCoord(lat.data, standard_name=None, long_name='latitude', var_name='latitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lat_coord
    
def make_lon_xstag_coordinate(filename):
    from iris import coords
    lon= loadwrfcube(filename, 'XLONG_U')
    lon_coord=coords.AuxCoord(lon.data, standard_name=None, long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lon_coord
    
def make_lat_xstag_coordinate(filename):
    from iris import coords
    lat= loadwrfcube(filename, 'XLAT_U')
    lat_coord=coords.AuxCoord(lat.data, standard_name=None, long_name='latitude', var_name='latitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lat_coord

def make_lon_ystag_coordinate(filename):
    from iris import coords
    lon= loadwrfcube(filename, 'XLONG_V')
    lon_coord=coords.AuxCoord(lon.data, standard_name=None, long_name='longitude', var_name='longitude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return lon_coord
    
def make_lat_ystag_coordinate(filename):
    from iris import coords
    lat= loadwrfcube(filename, 'XLAT_V')
    lat_coord=coords.AuxCoord(lat.data, standard_name=None, long_name='latitude', var_name='latidude', units='degrees', bounds=None, attributes=None, coord_system=None)
    return  lat_coord

