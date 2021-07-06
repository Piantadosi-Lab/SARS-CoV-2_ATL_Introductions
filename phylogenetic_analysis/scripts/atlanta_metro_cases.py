import pandas as pd
import datetime

dat_path = 'data/time_series_covid19_confirmed_US.csv'

dat = pd.read_csv(dat_path)

arc_counties = ['Cherokee', 'Clayton', 'Cobb', 'DeKalb', 'Douglas', 'Fayette', 'Fulton', 'Gwinnett', 'Henry', 'Rockdale']
arc_data = dat[(dat['Province_State']=='Georgia') & (dat['Admin2'].isin(arc_counties))]
arc_data.index = arc_data['Admin2']
arc_data = arc_data.iloc[:,11:]
arc_data.columns = pd.to_datetime(arc_data.columns)
arc_data = arc_data.loc[:,[i for i in arc_data.columns if i <= datetime.datetime.strptime('2020-03-31', '%Y-%m-%d')]]
arc_data.sum()['2020-03-31']

ga_data = dat[dat['Province_State'] == 'Georgia']
ga_data = ga_data.iloc[:,11:]
ga_data.columns = pd.to_datetime(ga_data.columns)
ga_data = ga_data.loc[:,[i for i in ga_data.columns if i <= datetime.datetime.strptime('2020-03-31', '%Y-%m-%d')]]


