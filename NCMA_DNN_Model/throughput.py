import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

file = 'D:\\DataSet\\单天线\\BPSK\\degree-SNR=8diff-A-B-C-D-XOR.mat'  # matlib文件位置
data = loadmat(file, mat_dtype=True)  # mat_dtype=True，保证了导入后变量的数据类型与原类型一致。
phy_raw_map = data['phy_raw_map']  # 导入后的data是一个字典，取出想要的变量字段即可。
# data_values = phy_raw_map[:, (0,1)]
# labels_values = np.delete(phy_raw_map, 0, axis=1)
# labels_values = np.delete(labels_values, 0, axis=1)
# 标签向量化
# x_values = np.asarray(data_values)
# # y_values = to_categorical(labels_values)
# y_values = np.asarray(labels_values)
r = len(phy_raw_map)
b = np.zeros((r,1))
phy_raw_map = np.c_[b,phy_raw_map]

# phy_raw_map = model_1.predict(x_values)
# phy_raw_map = np.asarray(phy_raw_map)

# b = np.zeros((r,1))
# phy_raw_map = np.c_[b,b,b,phy_raw_map]
# print(phy_raw_map.shape)
# print(phy_raw_map[:,9])
list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
for a in list:
    sum1 = 0
    sum2 = 0
    sum3 = 0
    sum4 = 0
    sum5 = 0
    sum6 = 0
    sum7 = 0
    sum8 = 0
    sum = 0
    upper_bound = 0
    i = 1
    b = a
    for i in range(r):
        if phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                phy_raw_map[i, 6] < a and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a:
            sum8 = sum8 + 1
        elif (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
              phy_raw_map[i, 6] > b and phy_raw_map[i, 7] > b and phy_raw_map[i, 8] > b and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] > b and phy_raw_map[i, 7] > b and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] > b and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] > b and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] > b and phy_raw_map[i, 8] > b and phy_raw_map[i, 9] < a):
            sum7 = sum7 + 1
        elif (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] > b and \
              phy_raw_map[i, 6] > b and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] > b and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] > b and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] > b and phy_raw_map[i, 4] > b and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] > b and phy_raw_map[i, 9] < a):
            sum6 = sum6 + 1
        elif (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
              phy_raw_map[i, 6] > b and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] > b and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] > b and phy_raw_map[i, 9] < a):
            sum5 = sum5 + 1
        elif phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                phy_raw_map[i, 6] < a and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] > b:
            sum4 = sum4 + 1
        elif (phy_raw_map[i, 3] > b and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and phy_raw_map[i, 6] < a and
              phy_raw_map[i, 7] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] > b and phy_raw_map[i, 5] < a and phy_raw_map[i, 6] < a and
                 phy_raw_map[i, 8] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] > b and phy_raw_map[i, 7] < a and
                 phy_raw_map[i, 8] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] > b and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] > b) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] > b and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] > b) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] > b and phy_raw_map[i, 9] > b):
            sum1 = sum1 + 1
        elif (phy_raw_map[i, 3] > b and phy_raw_map[i, 4] > b and phy_raw_map[i, 5] < a and \
              phy_raw_map[i, 7] < a and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] > b and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] > b and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] > b and phy_raw_map[i, 5] > b and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] > b and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] > b and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] > b and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] > b and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] > b and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] > b and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] > b and phy_raw_map[i, 5] < a and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] > b and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] > b and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] > b and phy_raw_map[i, 8] < a and phy_raw_map[i, 9] < a) or \
                (phy_raw_map[i, 3] < a and phy_raw_map[i, 4] < a and phy_raw_map[i, 5] > b and \
                 phy_raw_map[i, 6] < a and phy_raw_map[i, 7] < a and phy_raw_map[i, 8] > b and phy_raw_map[i, 9] < a):
            sum2 = sum2 + 1
        else:
            sum3 = sum3 + 1
    sum = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8
    upper_bound = (3 * sum3 + 2 * (sum2 + sum6 + sum7) + sum1 + sum4 + sum5) / r
    print(sum)
    print('upper_bound =', upper_bound)
    print('sum1 =', sum1, 'sum2 =', sum2, 'sum3 =', sum3, 'sum4 =', sum4, 'sum5 =', sum5, 'sum6 =', sum6, 'sum7 =',sum7, 'sum8 =', sum8)
