import matplotlib.pyplot as plt
import numpy as np
#from pylab import * #pylab 是 matplotlib 面向对象绘图库的一个接口

# XAy=[0,0]
# XBy=[0,0]
# XAB1y=[0,0]
# XAB1x=[0,0]
# XAx=[-1,1]
# XBx=[-1,1]
# XAB2x=[-2,2]
# XAB2y=[0,0]

XAy=[0,0]
XBy=[-1,1]
XAB1y=[-1,1]
XAB1x=[-1,1]
XAx=[-1,1]
XBx=[0,0]
XAB2x=[1,-1]
XAB2y=[-1,1]

plt.xlim(-2.5,2.5)# 设置横轴的上下限
plt.xticks(np.linspace(-2,2,5,endpoint=True))# 设置横轴记号, 是一个 numpy 数组，包含了从 −2 到 +2 等间隔的 5 个值
plt.ylim(-2.5,2.5)
plt.xticks(np.linspace(-2,2,5,endpoint=True))
plt.grid(linestyle=":")

# plt.scatter(XAx,XAy,s=150,c='w',linewidths=2,edgecolors='Black',marker="o",alpha=0.5,label="$X_A$")
# plt.scatter(XBx,XBy,s=150,c='w',linewidths=2,edgecolors='Black',marker="v",alpha=0.5,label="$X_B$")
# plt.scatter(XAB1x,XAB1y,s=150,c='w',linewidths=2,edgecolors='Black',marker="d",alpha=0.5,label="$X_A+X_B$")
# plt.scatter(XAB2x,XAB2y,s=150,c='w',linewidths=2,edgecolors='Black',marker="d",alpha=0.5)

plt.scatter(XAx,XAy,s=150,c='w',linewidths=2,edgecolors='Black',marker="o",alpha=0.5,label="$X_A$")
plt.scatter(XBx,XBy,s=150,c='w',linewidths=2,edgecolors='Black',marker="v",alpha=0.5,label="$X_Be^{j\pi/2}$")
plt.scatter(XAB1x,XAB1y,s=150,c='w',linewidths=2,edgecolors='Black',marker="d",alpha=0.5,label="$X_A+X_Be^{j\pi/2}$")
plt.scatter(XAB2x,XAB2y,s=150,c='w',linewidths=2,edgecolors='Black',marker="d",alpha=0.5)
plt.legend(loc="lower right")
plt.show()

