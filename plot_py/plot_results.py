import numpy as np
import matplotlib.pyplot as plt


N_train = 240
data_all = np.loadtxt("../cmake-build-debug/data/IPCLK_mean_data.txt")
data_pred_SARIMA = np.loadtxt("../cmake-build-debug/data/sarima_pred_test.txt")
data_pred_Linear_SARIMA = np.loadtxt("../cmake-build-debug/data/linear_sarima_pred_test.txt")
data_pred_test_linear = np.loadtxt("../cmake-build-debug/data/linear_pred_test.txt")
data_pred_train_linear = np.loadtxt("../cmake-build-debug/data/linear_pred_train.txt")
# data_pred_test_mean = np.loadtxt("../cmake-build-debug/data/mean_pred_test.txt")
# data_pred_train_mean = np.loadtxt("../cmake-build-debug/data/mean_pred_train.txt")

x_pred_test = np.arange(N_train,N_train+len(data_pred_SARIMA))
x_pred_train = np.arange(N_train)

plt.figure()
plt.plot(data_all,label='original mean data')
plt.plot(x_pred_test,data_pred_SARIMA,label='SARIMA predict')
plt.plot(x_pred_test,data_pred_Linear_SARIMA,label='Linear_SARIMA predict')
plt.plot(x_pred_train,data_pred_train_linear,label='linear fit train')
plt.plot(x_pred_test,data_pred_test_linear,label='linear predict test')
# plt.plot(x_pred_train,data_pred_train_mean,label='mean fit train')
# plt.plot(x_pred_test,data_pred_test_mean,label='mean predict test')
plt.legend()
plt.show()

