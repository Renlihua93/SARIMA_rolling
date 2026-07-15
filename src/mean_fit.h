#ifndef SKEW_CALC_H
#define SKEW_CALC_H

// 计算均值偏移
double calculate_mean_skew(double time_data[], int length);

// 计算偏置项 bias
double mean_fit(double mean_skew, double train_x[], double train_y[], int data_len);

// 线性预测
void mean_predict(double mean_skew, double bias,double *x, double*pred_y,int n);

#endif // SKEW_CALC_H