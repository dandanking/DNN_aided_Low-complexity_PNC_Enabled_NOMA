import numpy as np
import matplotlib.pyplot as plt

ax1 = plt.subplot(1, 1, 1)  # 第一行第一列图形
# ax2 = plt.subplot(1, 2, 2)  # 第一行第二列图形
list = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9']
width = 0.2
x = np.arange(len(list))  # the label locations
# dnnupper = [3.53747,
# 3.41532,
# 3.33646,
# 3.27463,
# 3.21912,
# 3.1646,
# 3.10568,
# 3.03353,
# 2.90216
# ]
difference1=[0.68605,
0.5639,
0.48504,
0.42321,
0.3677,
0.31318,
0.25426,
0.18211,
0.05074
]
difference2=[0.49423,
0.39208,
0.31322,
0.25139,
0.19588,
0.14136,
0.08244,
0.01029,
0.12108

]
difference3=[0.35397,
0.25182,
0.17296,
0.05782,
0.0011,
0.05562,
0.10113,

0.15997,
0.26134,

]
difference4=[0.26002,
0.13787,
0.05901,
0.00282,
0.05833,
0.11285,
0.17177,
0.24392,
0.37529
]
difference5=[0.16533,
0.04318,
0.03568,
0.09751,
0.15302,
0.20754,
0.26646,
0.33861,
0.46998

]
difference6=[0.08905,
0.0331,
0.11196,
0.17379,
0.2293,
0.28382,
0.34274,
0.41489,
0.54626
]
difference7=[
0.24035,
0.1182,
0.03934,
0.02249,
0.078,
0.13252,
0.19144,
0.26359,
0.39496
]
# list = [7,8,9,10]
# difference1 = [2.85142,
#                3.02324,
#                3.1635,
#                3.27745,
#                3.37214,
#                3.44842]
plt.sca(ax1)
plt.grid(linestyle=":")
plt.ylim(0,0.55)
# plt.plot([list[0], list[8]], [difference1[0], difference1[0]], '-', label='SNR=8(traditional)')
# plt.plot([list[0], list[8]], [difference1[1], difference1[1]], '-', label='SNR=9(traditional)')
# plt.plot([list[0], list[8]], [difference1[2], difference1[2]], '-', label='SNR=10(traditional)')
# plt.plot([list[0], list[8]], [difference1[3], difference1[3]], '-', label='SNR=11(traditional)')
# plt.plot([list[0], list[8]], [difference1[4], difference1[4]], '-', label='SNR=12(traditional)')
# plt.plot([list[0], list[8]], [difference1[5], difference1[5]], '-', label='SNR=13(traditional)')
# plt.plot(list, difference1[1], 'o-',label='SNR=8')
# plt.plot(list, difference1[2], 'o-',label='SNR=9')
# plt.plot(list, difference1[3], 'o-',label='SNR=10')
# plt.plot(list, dnnupper, 'o-', label='SNR=10(DNN-Aided)')
#plt.plot(list, difference1, 'o-', label='SNR=8')
plt.bar(x , difference7,width)#, label='SNR=9'
#plt.bar(x, difference3, width)#, label='SNR=10'
#plt.bar(x, difference4, width)#, label='SNR=11'
#plt.bar(x, difference5, width)#, label='SNR=12'
#plt.plot(list, difference6, 'o-', label='SNR=13')
#plt.title('Throughout(packets per slot)')
plt.xticks(x,list)
plt.xlabel('$\\alpha$')
plt.ylabel('Throughput difference (packets per slot)')
#plt.legend()
plt.show()
