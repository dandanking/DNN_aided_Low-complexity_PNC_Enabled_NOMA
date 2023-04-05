import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from pandas import DataFrame

file = 'D:\\DataSet\\单天线\\BPSK\\degree-diffSNR-A-B-C-D-XOR.mat'  # matlib文件位置
data = loadmat(file, mat_dtype=True)  # mat_dtype=True，保证了导入后变量的数据类型与原类型一致。
phy_raw_map = data['phy_raw_map']  # 导入后的data是一个字典，取出想要的变量字段即可。
#data_values = phy_raw_map[:, (0,1)]
#labels_values = np.delete(phy_raw_map, 0, axis=1)
#labels_values = np.delete(labels_values, 0, axis=1)
# 标签向量化
#x_values = np.asarray(data_values)
# y_values = to_categorical(labels_values)
#y_values = np.asarray(labels_values)
########################################################################################
#phy_raw_map = model_1.predict(x_values)
phy_raw_map = np.asarray(phy_raw_map)
r = len(phy_raw_map)
b = np.zeros((r,1))
phy_raw_map = np.c_[b,phy_raw_map]
print(phy_raw_map.shape)
print(phy_raw_map[:,18])
list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
a=0.1
sum1 = 0
sum2 = 0
sum3 = 0
sum4 = 0
sum = 0
upper_bound = 0
i = 1
for i in range(r):
        if phy_raw_map[i, 8] >a and phy_raw_map[i, 14] >a or phy_raw_map[i, 10]>a and phy_raw_map[i, 16]>a or \
                phy_raw_map[i, 12]>a and phy_raw_map[i, 17]>a or phy_raw_map[i, 15]>a and phy_raw_map[i, 18]>a or \
                phy_raw_map[i, 4] > a and phy_raw_map[i, 9] > a or phy_raw_map[i, 5]>a and phy_raw_map[i, 11]>a or \
                phy_raw_map[i, 7] > a and phy_raw_map[i, 13] > a or phy_raw_map[i, 4]>a and phy_raw_map[i, 5]>a and phy_raw_map[i, 14]>a or \
                phy_raw_map[i, 4]>a and phy_raw_map[i, 7]>a and phy_raw_map[i, 16]>a or phy_raw_map[i, 5]>a and phy_raw_map[i, 7]>a and phy_raw_map[i, 17]>a or \
                phy_raw_map[i, 7]>a and phy_raw_map[i, 8]>a and phy_raw_map[i, 18]>a or phy_raw_map[i, 5]>a and phy_raw_map[i, 10]>a and phy_raw_map[i, 18]>a or \
                phy_raw_map[i, 4]>a and phy_raw_map[i, 12]>a and phy_raw_map[i, 18]>a or phy_raw_map[i, 4]>a and phy_raw_map[i, 5]>a and phy_raw_map[i, 7]>a and phy_raw_map[i, 18] or \
                phy_raw_map[i, 6]>a:
            phy_raw_map[i, 6]=1
        if phy_raw_map[i, 7]>a or phy_raw_map[i, 4] >a and phy_raw_map[i, 10] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 12] >a or phy_raw_map[i, 6] >a and phy_raw_map[i, 13] >a or \
                phy_raw_map[i, 8] > a and phy_raw_map[i, 15] >a or phy_raw_map[i, 9] > a and phy_raw_map[i, 16] > a or phy_raw_map[i, 11] >a and phy_raw_map[i, 17]>a or \
                phy_raw_map[i, 14] > a and phy_raw_map[i, 18]>a or phy_raw_map[i, 4] > a and phy_raw_map[i, 5]>a and phy_raw_map[i, 15] >a or \
                phy_raw_map[i, 4] >a and phy_raw_map[i, 6]>a and phy_raw_map[i, 16] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 6]>a and phy_raw_map[i, 17] >a or \
                phy_raw_map[i, 6] > a and phy_raw_map[i, 8] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 9]>a and phy_raw_map[i, 18] >a or \
                phy_raw_map[i, 4] > a and phy_raw_map[i, 11] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 5]>a and phy_raw_map[i, 6] >a and phy_raw_map[i, 18]>a:
            phy_raw_map[i, 7]=1
        if phy_raw_map[i, 4]>a or phy_raw_map[i, 5] >a and phy_raw_map[i, 8] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 9] >a or phy_raw_map[i, 7] >a and phy_raw_map[i, 10] >a or \
                phy_raw_map[i, 11] > a and phy_raw_map[i, 14] >a or phy_raw_map[i, 12] > a and phy_raw_map[i, 15] > a or phy_raw_map[i, 13] >a and phy_raw_map[i, 16]>a or \
                phy_raw_map[i, 17] > a and phy_raw_map[i, 18]>a or phy_raw_map[i, 5] > a and phy_raw_map[i, 6]>a and phy_raw_map[i, 14] >a or \
                phy_raw_map[i, 5] >a and phy_raw_map[i, 7]>a and phy_raw_map[i, 15] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 7]>a and phy_raw_map[i, 16] >a or \
                phy_raw_map[i, 7] > a and phy_raw_map[i, 11] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 12]>a and phy_raw_map[i, 18] >a or \
                phy_raw_map[i, 5] > a and phy_raw_map[i, 13] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 6]>a and phy_raw_map[i, 7] >a and phy_raw_map[i, 18]>a:
            phy_raw_map[i, 4]=1
        if phy_raw_map[i, 5]>a or phy_raw_map[i, 4] >a and phy_raw_map[i, 8] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 11] >a or phy_raw_map[i, 7] >a and phy_raw_map[i, 12] >a or \
                phy_raw_map[i, 9] > a and phy_raw_map[i, 14] >a or phy_raw_map[i, 10] > a and phy_raw_map[i, 15] > a or phy_raw_map[i, 13] >a and phy_raw_map[i, 17]>a or \
                phy_raw_map[i, 16] > a and phy_raw_map[i, 18]>a or phy_raw_map[i, 4] > a and phy_raw_map[i, 6]>a and phy_raw_map[i, 14] >a or \
                phy_raw_map[i, 4] >a and phy_raw_map[i, 7]>a and phy_raw_map[i, 15] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 7]>a and phy_raw_map[i, 17] >a or \
                phy_raw_map[i, 7] > a and phy_raw_map[i, 9] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 10]>a and phy_raw_map[i, 18] >a or \
                phy_raw_map[i, 4] > a and phy_raw_map[i, 13] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 6]>a and phy_raw_map[i, 7] >a and phy_raw_map[i, 18]>a:
            phy_raw_map[i, 5]=1
        if phy_raw_map[i, 8]>a or phy_raw_map[i, 4] >a and phy_raw_map[i, 5] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 14] > a or phy_raw_map[i, 7] >a and phy_raw_map[i, 15] > a or \
                phy_raw_map[i, 13] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 7] > a and phy_raw_map[i, 18] > a:
            phy_raw_map[i, 8]=1
        if phy_raw_map[i, 9]>a or phy_raw_map[i, 4] >a and phy_raw_map[i, 6] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 14] > a or phy_raw_map[i, 7] >a and phy_raw_map[i, 16] > a or \
                phy_raw_map[i, 12] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 7] > a and phy_raw_map[i, 18] > a:
            phy_raw_map[i, 9]=1
        if phy_raw_map[i, 10]>a or phy_raw_map[i, 4] >a and phy_raw_map[i, 7] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 15] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 16] > a or \
                phy_raw_map[i, 11] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 6] > a and phy_raw_map[i, 18] > a:
            phy_raw_map[i, 10]=1
        if phy_raw_map[i, 11]>a or phy_raw_map[i, 5] >a and phy_raw_map[i, 6] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 14] > a or phy_raw_map[i, 7] >a and phy_raw_map[i, 17] > a or \
                phy_raw_map[i, 10] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 7] > a and phy_raw_map[i, 18] > a:
            phy_raw_map[i, 11]=1
        if phy_raw_map[i, 12]>a or phy_raw_map[i, 5] >a and phy_raw_map[i, 7] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 15] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 17] > a or \
                phy_raw_map[i, 9] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 6] > a and phy_raw_map[i, 18] > a:
            phy_raw_map[i, 12]=1
        if phy_raw_map[i, 13]>a or phy_raw_map[i, 6] >a and phy_raw_map[i, 7] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 16] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 17] > a or \
                phy_raw_map[i, 8] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 5] > a and phy_raw_map[i, 18] > a:
            phy_raw_map[i, 13]=1
        if phy_raw_map[i, 14]>a or phy_raw_map[i, 6] >a and phy_raw_map[i, 8] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 9] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 11] > a or \
                phy_raw_map[i, 7] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 13] >a and phy_raw_map[i, 15] > a or phy_raw_map[i, 12] >a and phy_raw_map[i, 16] > a or \
                phy_raw_map[i, 10] > a and phy_raw_map[i, 17] > a:
            phy_raw_map[i, 14]=1
        if phy_raw_map[i, 15]>a or phy_raw_map[i, 7] >a and phy_raw_map[i, 8] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 10] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 12] > a or \
                phy_raw_map[i, 6] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 13] >a and phy_raw_map[i, 14] > a or phy_raw_map[i, 11] >a and phy_raw_map[i, 16] > a or \
                phy_raw_map[i, 9] > a and phy_raw_map[i, 17] > a:
            phy_raw_map[i, 15]=1
        if phy_raw_map[i, 16]>a or phy_raw_map[i, 7] >a and phy_raw_map[i, 9] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 10] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 13] > a or \
                phy_raw_map[i, 5] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 12] >a and phy_raw_map[i, 14] > a or phy_raw_map[i, 11] >a and phy_raw_map[i, 15] > a or \
                phy_raw_map[i, 8] > a and phy_raw_map[i, 17] > a:
            phy_raw_map[i, 16]=1
        if phy_raw_map[i, 17]>a or phy_raw_map[i, 7] >a and phy_raw_map[i, 11] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 12] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 13] > a or \
                phy_raw_map[i, 4] > a and phy_raw_map[i, 18] > a or phy_raw_map[i, 10] >a and phy_raw_map[i, 14] > a or phy_raw_map[i, 9] >a and phy_raw_map[i, 15] > a or \
                phy_raw_map[i, 8] > a and phy_raw_map[i, 16] > a:
            phy_raw_map[i, 17]=1
        if phy_raw_map[i, 18]>a or phy_raw_map[i, 8] >a and phy_raw_map[i, 13] > a or phy_raw_map[i, 9] >a and phy_raw_map[i, 12] > a or phy_raw_map[i, 10] >a and phy_raw_map[i, 11] > a or \
                phy_raw_map[i, 7] > a and phy_raw_map[i, 14] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 15] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 16] > a or \
                phy_raw_map[i, 4] >a and phy_raw_map[i, 17] > a or phy_raw_map[i, 6] >a and phy_raw_map[i, 7] > a and phy_raw_map[i, 8] > a or \
                phy_raw_map[i, 5] >a and phy_raw_map[i, 7] > a and phy_raw_map[i, 9] > a or phy_raw_map[i, 5] >a and phy_raw_map[i, 6] > a and phy_raw_map[i, 10] > a or \
                phy_raw_map[i, 4] > a and phy_raw_map[i, 7] > a and phy_raw_map[i, 11] > a or phy_raw_map[i, 4] >a and phy_raw_map[i, 6] > a and phy_raw_map[i, 12] > a or \
                phy_raw_map[i, 4] > a and phy_raw_map[i, 5] > a and phy_raw_map[i, 13] > a:
            phy_raw_map[i, 18]=1

for i in range(r):
        if phy_raw_map[i, 4]>a and phy_raw_map[i,5]>a and phy_raw_map[i,6]>a and phy_raw_map[i,7]>a:
            sum1=sum1+1
        elif phy_raw_map[i, 4]>a and phy_raw_map[i,5]>a and phy_raw_map[i,6]>a or phy_raw_map[i, 4]>a and phy_raw_map[i,5]>a and phy_raw_map[i,7]>a or \
                phy_raw_map[i, 4]>a and phy_raw_map[i,6]>a and phy_raw_map[i,7]>a or phy_raw_map[i, 5]>a and phy_raw_map[i,6]>a and phy_raw_map[i,7]>a:
            sum2=sum2+1
        elif phy_raw_map[i,4]>a and phy_raw_map[i,5]>a or phy_raw_map[i,4]>a and phy_raw_map[i,6]>a or phy_raw_map[i,4]>a and phy_raw_map[i,7]>a or \
                phy_raw_map[i,5]>a and phy_raw_map[i,6]>a or phy_raw_map[i,5]>a and phy_raw_map[i,7]>a or phy_raw_map[i,6]>a and phy_raw_map[i,7]>a:
            sum3=sum3+1
        else:
            sum4=sum4+1
#df=DataFrame(phy_raw_map)
#df.to_csv("D:\DataSet//df.csv")
sum = sum1 + sum2 + sum3 + sum4
upper_bound = (4 * sum1 + 3 * sum2 + 3*sum3 + 2*sum4 ) / r
print(sum)
print('upper_bound =', upper_bound)
print('sum1 =', sum1, 'sum2 =', sum2, 'sum3 =', sum3, 'sum4 =', sum4)
