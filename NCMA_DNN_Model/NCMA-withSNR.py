import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
from tensorflow.keras import layers
from scipy.io import loadmat
from tensorflow.keras.utils import to_categorical

# 准备数据阶段
file = 'D:\\DataSet\degree-SNR-A-B-XOR.mat'  # matlib文件位置
data = loadmat(file, mat_dtype=True)  # mat_dtype=True，保证了导入后变量的数据类型与原类型一致。
phy_raw_map = data['phy_raw_map']  # 导入后的data是一个字典，取出想要的变量字段即可。
data_values = phy_raw_map[:, 0]
labels_values = np.delete(phy_raw_map, 0, axis=1)
labels_values = np.delete(labels_values, 0, axis=1)

# 标签向量化
x_values = np.asarray(data_values)
# y_values = to_categorical(labels_values)
y_values = np.asarray(labels_values)

# 划分数据集
SAMPLES = 100000
TRAIN_SPLIT = int(0.6 * SAMPLES)
TEST_SPLIT = int(0.2 * SAMPLES + TRAIN_SPLIT)
x_train, x_validate, x_test = np.split(x_values, [TRAIN_SPLIT, TEST_SPLIT])
y_train, y_validate, y_test = np.split(y_values, [TRAIN_SPLIT, TEST_SPLIT])
assert (x_train.size + x_validate.size + x_test.size) == SAMPLES
print(x_train.shape, x_validate.shape, x_test.shape)
print(y_train.shape, y_validate.shape, y_test.shape)

# 模型
model_1 = tf.keras.Sequential()
model_1.add(layers.Dense(50, activation='relu', input_shape=(1,)))
model_1.add(layers.Dense(50, activation='relu'))
model_1.add(layers.Dense(50, activation='relu'))
model_1.add(layers.Dense(50, activation='relu'))
model_1.add(layers.Dense(3, activation="sigmoid"))
model_1.compile(loss='binary_crossentropy', optimizer='rmsprop',  metrics=['mae'])
model_1.summary()

# 训练
history_1 = model_1.fit(x_train, y_train, epochs=500, batch_size=64, validation_data=(x_validate, y_validate))

# 绘制历史数据
loss = history_1.history['loss']
val_loss = history_1.history['val_loss']
mae = history_1.history['mae']
val_mae = history_1.history['val_mae']
epochs = range(1, len(loss) + 1)
predictions = model_1.predict(x_test)

ax1 = plt.subplot(2, 2, 1)  # 第一行第一列图形
ax2 = plt.subplot(2, 2, 2)  # 第一行第二列图形
ax3 = plt.subplot(2, 1, 2)  # 第二行
# ax4 = plt.subplot(3, 1, 3)  # 第三行
# ax5 = plt.subplot(2, 1, 2)  # 第四行

SKIP = 20
plt.sca(ax1)
plt.plot(epochs[SKIP:], loss[SKIP:], 'g.', label='Training loss')
plt.plot(epochs[SKIP:], val_loss[SKIP:], 'b.', label='Validation loss')
plt.title('Training and validation loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()

plt.sca(ax2)
plt.plot(epochs[SKIP:], mae[SKIP:], 'g.', label='Training mae')
plt.plot(epochs[SKIP:], val_mae[SKIP:], 'b.', label='Validation mae')
plt.title('Training and validation mae')
plt.xlabel('Epochs')
plt.ylabel('mae')
plt.legend()

plt.sca(ax3)
plt.plot(x_test, predictions[:, 0], 'g.', label='predictions A')
plt.plot(x_test, y_test[:, 0], 'r.', label="real A")
plt.title('real A and prediction A')
plt.xlabel('phase offsets')
plt.ylabel('probability')
plt.legend()
plt.show()

# plt.sca(ax4)
plt.plot(x_test, predictions[:, 1], 'g.', label='predictions B')
plt.plot(x_test, y_test[:, 1], 'r.', label="real B")
plt.title('real B and prediction B')
plt.xlabel('phase offsets')
plt.ylabel('probability')
plt.legend()


# plt.sca(ax5)
# plt.plot(x_test, predictions[:, 2], 'g.', label='predictions PNC')
# plt.plot(x_test, y_test[:, 2], 'r.', label="real PNC")
# plt.title('real PNC and prediction PNC')
# plt.xlabel('phase offsets')
# plt.ylabel('probability')
# plt.legend()
plt.show()

test1 = np.arange(0, 10)
test2 = np.arange(175, 185)
test3 = np.arange(350, 360)
predictions1 = model_1.predict(test1)
predictions2 = model_1.predict(test2)
predictions3 = model_1.predict(test3)
print(predictions1)
print(predictions2)
print(predictions3)
