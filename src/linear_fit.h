//
// Created by root on 2026/5/14.
//

#ifndef SARIAM_C_LINEAR_FIT_H
#define SARIAM_C_LINEAR_FIT_H

void linear_fit(
    const double *x,
    const double *y,
    int n,
    double *a,   // 斜率
    double *b    // 截距
);

void linear_predict(
    double k,   // 斜率
    double b,    // 截距
    double *x,
    double *pred_y,
    int n
);
#endif //SARIAM_C_LINEAR_FIT_H