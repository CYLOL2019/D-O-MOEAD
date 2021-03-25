import tushare as ts
import datetime
import numpy as np
import pandas as pd
import time
import collections
"""
"""
# Get the return of the stocks market A, SH and SZ, from 20181130 to 20201130.

# Load token
T_file = 'TOKEN.txt'
with open(T_file, 'r') as fin:
    token = fin.readline()
ts.set_token(token)
pro = ts.pro_api()
startime = '20090101'
endtime = '20201211'
# Get all_ts_codes
all_ts_code = pro.query('stock_basic', exchange='', list_status='L',
                        fields='ts_code,symbol,name,area,industry,list_date')
# Number of assets, day interval, close data
NoA = all_ts_code.shape[0]
# Number of workday
workday_df = pro.trade_cal(exchange='SSE', is_open='1',
                           start_date=startime, end_date=endtime,
                           fields='cal_date')
maxdays = workday_df.size
close_df = pd.DataFrame(np.full([maxdays, NoA], np.nan))
# print(all_ts_code['ts_code'].values)
close_df.columns = all_ts_code['ts_code'].values
# print(close_df)
close_df.index = workday_df['cal_date'].values[::-1]
# print(close_df)
print(maxdays)
for _ in all_ts_code['ts_code'].values:
    # print(all_ts_code)
    # print(_)
    # print(_)
    time.sleep(0.33)  # 200 times every minute
    data = pro.daily_basic(ts_code=_,
                           start_date=startime, end_date=endtime,
                           fields='close')
    print(data)
    if data.size == maxdays:
        data.index = workday_df['cal_date'].values[::-1]
        print(data)
        close_df[_] = data
    print(close_df)


# Truncation
values_set = close_df.iloc[0].values
index = []
all_ts_code = all_ts_code['ts_code'].values
for _ in range(NoA):
    if np.isnan(values_set[_]):
        index = np.union1d(index, all_ts_code[_])
for _ in index:
    close_df.drop(_, axis=1, inplace=True)
dirname = startime + ',' + endtime
outputfile = dirname + '/'  "data_tushare_" + endtime + ".csv"
close_df.to_csv(outputfile, index=True, sep=',')




