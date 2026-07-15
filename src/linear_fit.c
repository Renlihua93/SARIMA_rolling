//
// Created by root on 2026/5/14.
//

#include <stdio.h>
#include <stdlib.h>
#include "linear_fit.h"

void linear_fit(
    const double *x,
    const double *y,
    int n,
    double *a,   // 斜率
    double *b    // 截距
) {
    double sum_x = 0.0, sum_y = 0.0;
    double sum_xx = 0.0, sum_xy = 0.0;

    for (int i = 0; i < n; i++) {
        sum_x  += x[i];
        sum_y  += y[i];
        sum_xx += x[i] * x[i];
        sum_xy += x[i] * y[i];
    }

    double denom = n * sum_xx - sum_x * sum_x;
    if (denom == 0.0) {
        *a = 0.0;
        *b = 0.0;
        return;
    }

    *a = (n * sum_xy - sum_x * sum_y) / denom;
    *b = (sum_y - (*a) * sum_x) / n;
}

void linear_predict(
    double k,   // 斜率
    double b,    // 截距
    double *x,
    double *pred_y,
    int n
)
{
    for (int i = 0; i < n; i++)
    {
        pred_y[i] = k * x[i]+b;
    }

}