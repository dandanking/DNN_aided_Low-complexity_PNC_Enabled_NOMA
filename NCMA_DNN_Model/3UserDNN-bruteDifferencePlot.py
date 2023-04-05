import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
ax1 = plt.subplot(1, 1, 1)  # 第一行第一列图形
#ax2 = plt.subplot(1, 2, 2)  # 第一行第二列图形
list = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9']
width = 0.17
x = np.arange(len(list))  # the label locations
# dnnupper=[2.35016,
# 2.27118,
# 2.2199,
# 2.17486,
# 2.13628,
# 2.0975,
# 2.05536,
# 2.00458,
# 1.9282
# ]
difference1=[0.37181,
0.29283,
0.24155,
0.19651,
0.15793,
0.11915,
0.07701,
0.02623,
0.05015
]
difference2=[0.22615,
0.14717,
0.09589,
0.05085,
0.01227,
0.02651,
0.06865,
0.11943,
0.19581,
]
difference3=[0.09962,
0.02064,
0.03064,
0.07568,
0.11426,
0.15304,
0.19518,
0.24596,
0.32234,
]
difference4=[0.00322,
0.0822,
0.13348,
0.17852,
0.2171,
0.25588,
0.29802,
0.3488,
0.42518,
]

differencediff=[0.30278,
0.2238,
0.17252,
0.12748,
0.0889,
0.05012,
0.00798,
0.0428,
0.11918
]

Reduced_complexity=[
    0.384,
    0.343,
    0.299,
    0.273
]
SNR=[
    7,8,9,10
]

#list = [7,8,9,10]
# difference1=[1.97835,
#              2.12401,
#              2.25054,
#              2.35338]
plt.sca(ax1)

#plt.xlim(15,19)# 设置横轴的上下限
plt.xticks(np.linspace(7,10,4,endpoint=True))# 设置横轴记号, 是一个 numpy 数组，包含了从 7 到 10 等间隔的 5 个值
plt.ylim(0,0.5)
plt.grid(linestyle=":")
#plt.plot([list[0],list[8]], [difference1[0],difference1[0]],'-',label='SNR=7(traditional)')
#plt.plot([list[0],list[8]], [difference1[1],difference1[1]],'-',label='SNR=8(traditional)')
#plt.plot([list[0],list[8]], [difference1[2],difference1[2]],'-',label='SNR=9(traditional)')
#plt.plot([list[0],list[8]], [difference1[3],difference1[3]],'-',label='SNR=10(traditional)')
#plt.bar(x, difference1, width) #,label='SNR=7'
#plt.bar(x, difference2, width) #,label='SNR=8'
#plt.bar(x , difference3, width) #,label='SNR=9'
#plt.bar(x, difference4, width) # ,label='SNR=10'

plt.plot(SNR, Reduced_complexity)
# plt.bar(x + 2.5*width, differencediff, width,label='SNRA,SNRB,SNRC are different')

#plt.xticks(x,list)
# plt.xlabel('$\\alpha$')
# plt.ylabel("Throughput difference (packets per slot)")
plt.xlabel("SNR (dB)")
plt.ylabel('Reduced complexity (percentage)')

def to_percent(temp, position):

  return "%1.0f"%(100*temp) + "%"

plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))

#plt.plot(list, dnnupper, 'o-', label='SNR=8(DNN-Aided)')
#plt.plot(list, difference3, 'o-', label='SNR=9')
#plt.plot(list, difference4, 'o-', label='SNR=10')
# plt.title('Throughout(packets per slot)')

#plt.legend()
plt.show()