import tushare as ts
import datetime
import numpy as np
import pandas as pd
import math
import time
import collections

"""
"""


# get portfolio data of every asset
# startime= '20181230'
# traintime = '20200227'
# endtime = '20201130'
startime= '20181211'
traintime = '20200201'
endtime = '20201210'
# Get all_ts_codes
dirname = startime + ',' + endtime
file = dirname + '/' + 'data_tushare_' + endtime + '.csv'
data = pd.read_csv(file, index_col=0)
all_ts_code = data.columns.values
NoA = all_ts_code.shape[0]
# Truncation train data
data = data.loc[traintime:startime]
print(data)
# Create u and cov csv
data_u = np.zeros(NoA)
data_u = (data.mean(axis=0) - data.iloc[-1, :]) / data.iloc[-1, :]
data_cov = (data-data.iloc[-1, :]).cov()
data_cov.index = all_ts_code
file = dirname + '/' + 'portAC_u_' + traintime + '.csv'
with open(file, 'w') as fout:
    fout.write(str(data.shape[1])+'\n')
    for _ in data_u:
        fout.write(str(_)+'\n')
outputfile = dirname + '/' + "portAC_cov_" + traintime + ".csv"
data_cov.to_csv(outputfile, index=True, sep=',')
# Create vailable_ts_code
# all_ts_code = data.columns.values
file = dirname + '/' + "available_ts_code_" + traintime + ".csv"
with open(file, 'w') as fout:
    fout.write('ts_code'+'\n')
    for _ in all_ts_code:
        fout.write(str(_)+'\n')


# portfolioCA = 'portCA.txt'
# with open(portfolioCA, 'w') as fout:
#     fout.write(data.shape[0]+'\n')





