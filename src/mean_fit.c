//
// Created by root on 2026/5/28.
//

#include <stdio.h>
#include <stdlib.h>


// time_data: 输入数组
// length: 数组长度（必须>=2，否则无法计算差值）
double calculate_mean_skew(double time_data[], int length) {
    // 数组长度不足，返回0（避免越界）
    if (length < 2) {
        return 0.0;
    }

    double sum = 0.0;
    // 计算相邻元素差值并求和
    for (int i = 0; i < length - 1; i++) {
        sum += (time_data[i + 1] - time_data[i])/0.02;
    }

    // 计算差值的平均值
    double mean_skew = sum / (length - 1);
    return mean_skew;
}

// 计算偏置项
double mean_fit(double mean_skew, double train_x[], double train_y[], int data_len) {
    // 计算train_x和train_y的平均值
    double sum_x = 0.0, sum_y = 0.0;
    for (int i = 0; i < data_len; i++) {
        sum_x += train_x[i];
        sum_y += train_y[i];
    }

    double mean_x = sum_x / data_len;
    double mean_y = sum_y / data_len;

    // 计算偏置 bias = mean(y) - k * mean(x)
    double bias = mean_y - mean_skew * mean_x;
    return bias;
}

// 预测函数
void mean_predict(double mean_skew, double bias,double *x, double*pred_y,int n) {
    for (int i = 0; i < n; i++)
    {
        pred_y[i] = mean_skew * x[i]+bias;
    }

}
