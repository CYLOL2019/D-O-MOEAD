# ICDV-MOEAD
At first, that is very kind of the owners of the financial data platform Tushare (https://tushare.pro). They provide free access for us.

datadownload.py introduces a very simple way to download all the data about close prices of every assets in Shanghai and Shenzhen stock markets from 01/01/2009 to 11/12/2020.
getport_ICDV.py introduces a program to convert 'data_tushare_20201210.csv' into 'portAC_u_20200201.csv', 'portAC_cov_20200201.csv' and 'available_ts_code_20200201.csv'.

The folder '20181211,2020201' contains the close prices of every assets from 11/12/2018 to 10/12/2020, 'data_tushare_20201210.csv'.
Moreover, the data from 11/12/2018 to 01/02/2020 is extracted to construct the expected return, 'portAC_u_20200201.csv', and risk.
Since the  size of 'portAC_cov_20200201.csv' reaches about 150MB, we can not upload it. One can form it with the details provided in the paper.
'available_ts_code_20200201.csv' contains all the code of the assets.
'000001.SH_20201210.csv' and '399001.SZ_20201210.csv' include the indices of Shanghai and Shenzhen stock markets respectively.
