import numpy as np
import pandas as pd
import pytest
import xarray as xr
from dssatutils.weather_gridded_common import convert_units, extract_netcdf_series

@pytest.mark.parametrize('value,unit', [(20,'MJ m-2 day-1'),(20000,'kJ m-2 day-1'),(20000000,'J m-2'),(20/.0864,'W m-2')])
def test_equivalent_daily_energy(value, unit):
    assert convert_units([value],unit,'srad')[0] == pytest.approx(20)

@pytest.mark.parametrize('kind,unit', [('srad','unknown'),('temp',''),('vp','unknown')])
def test_ambiguous_units_fail(kind,unit):
    with pytest.raises(ValueError, match='units'): convert_units([20],unit,kind)

def test_pressure_units_and_wind_height():
    assert convert_units([1200],'Pa','vp') == pytest.approx(convert_units([12],'hPa','vp'))
    assert convert_units([1.2],'kPa','vp') == pytest.approx(convert_units([12],'hPa','vp'))
    assert convert_units([3.6],'km/h','wind',wind_height_m=2)[0] == 1
    with pytest.raises(ValueError,match='height_m'): convert_units([1],'m/s','wind')

def dataset(tmp_path,times,axis_units='degrees_north'):
    ds=xr.Dataset({'tasmax':(('time','lat','lon'),np.full((len(times),2,2),300.),{'units':'K'})},
        coords={'time':times,'lat':('lat',[0.,1.],{'units':axis_units}),'lon':('lon',[0.,1.],{'units':'degrees_east'})})
    f=tmp_path/'tmax.nc';ds.to_netcdf(f);return f

def test_domain_is_not_extended_beyond_grid_centres(tmp_path):
    f=dataset(tmp_path,pd.date_range('2020-01-01',periods=2))
    with pytest.raises(ValueError,match='outside grid'):
        extract_netcdf_series(f,['tasmax'],['p'],[1.1],[.5],2020,2020,'temp')

def test_projected_axes_rejected(tmp_path):
    f=dataset(tmp_path,pd.date_range('2020-01-01',periods=2),'m')
    with pytest.raises(ValueError,match='geographic'):
        extract_netcdf_series(f,['tasmax'],['p'],[.5],[.5],2020,2020,'temp')

def test_subdaily_never_overwrites(tmp_path):
    f=dataset(tmp_path,pd.date_range('2020-01-01',periods=2,freq='6h'))
    with pytest.raises(ValueError,match='Subdaily'):
        extract_netcdf_series(f,['tasmax'],['p'],[.5],[.5],2020,2020,'temp')
