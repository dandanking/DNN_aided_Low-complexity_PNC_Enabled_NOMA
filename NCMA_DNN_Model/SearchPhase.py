import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

file = 'D:\\DataSet\degree-SNR=9-A-B-C.mat'  # matlib文件位置
data = loadmat(file, mat_dtype=True)  # mat_dtype=True，保证了导入后变量的数据类型与原类型一致。
phy_raw_map = data['phy_raw_map']  # 导入后的data是一个字典，取出想要的变量字段即可。
phy_raw_map = np.asarray(phy_raw_map)
#data_values = phy_raw_map[:, (0,1)]
#labels_values = np.delete(phy_raw_map, 0, axis=1)
#labels_values = np.delete(labels_values, 0, axis=1)
x=list()
y=list()
r = len(phy_raw_map)
for i in range(r):
    if phy_raw_map[i, 2]==0 and phy_raw_map[i, 3]==0 and phy_raw_map[i, 4]==0 and phy_raw_map[i, 5] and phy_raw_map[i, 6] and phy_raw_map[i, 7] and phy_raw_map[i, 8]:
       x.append(phy_raw_map[i, 0])
       y.append(phy_raw_map[i, 1])
       print(phy_raw_map[i, 0])
       print(phy_raw_map[i, 1])
plt.plot(x, y, 'r.')
plt.xlabel('Phase Offset B')
plt.ylabel('Phase Offset C')
#plt.legend()
plt.show()