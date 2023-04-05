import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from pandas import DataFrame

file = 'D:\\DataSet\\单天线\\QPSK\\2User-SNR=19.mat'  # matlib文件位置
data = loadmat(file, mat_dtype=True)  # mat_dtype=True，保证了导入后变量的数据类型与原类型一致。
phy_raw_map = data['phy_raw_map']  # 导入后的data是一个字典，取出想要的变量字段即可。
r = len(phy_raw_map)
phy_raw_map = np.asarray(phy_raw_map)
print(phy_raw_map.shape)
sum=0
for i in range(r):
    if phy_raw_map[i,0]==1:
        sum=sum+1
print(sum/5000)