// SPDX-License-Identifier: BSD-3-Clause
/*
 * sarimatest.c - SARIMA 滚动窗口示例程序
 *
 * 功能:
 * 1. 读取 data/IPCLK.txt
 * 2. 提取鉴相值
 * 3. 转换为相对鉴相值，单位 ns
 * 4. 做小时平均，降采样为 1 小时一个点
 * 5. 使用滚动窗口：
 *      起点: 第0天 到 第12天
 *      步长: 0.25天 = 6小时
 *      训练窗口: 10天
 *      预测窗口: 从训练结束预测到数据末尾
 * 6. 对每个窗口计算:
 *      Linear
 *      SARIMA
 *      Linear_SARIMA
 * 7. 保存预测结果和误差结果
 *
 * 须用 CMake 目标 sarimatest 链接 sarimalib 构建。
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "sarima.h"
#include "load_matrix.h"
#include "linear_fit.h"
#include "mean_fit.h"


/*
 * 计算误差指标
 *
 * 误差定义:
 *      error = actual - pred
 *
 * 输出:
 *      mae              平均绝对误差
 *      rmse             均方根误差
 *      max_abs_error    最大绝对误差
 *      exceed_100_count 绝对误差超过100ns的点数
 */
static void calc_error_metrics(
    const double *actual,
    const double *pred,
    int n,
    double *mae,
    double *rmse,
    double *max_abs_error,
    int *exceed_100_count
) {
    double sum_abs = 0.0;
    double sum_sq = 0.0;
    double max_abs = 0.0;
    int exceed_count = 0;

    for (int i = 0; i < n; i++) {
        double err = actual[i] - pred[i];
        double abs_err = fabs(err);

        sum_abs += abs_err;
        sum_sq += err * err;

        if (abs_err > max_abs) {
            max_abs = abs_err;
        }

        if (abs_err > 100.0) {
            exceed_count++;
        }
    }

    if (n > 0) {
        *mae = sum_abs / n;
        *rmse = sqrt(sum_sq / n);
    } else {
        *mae = 0.0;
        *rmse = 0.0;
    }

    *max_abs_error = max_abs;
    *exceed_100_count = exceed_count;
}


int main(void) {
    int i;

    /*
     * =========================
     * 1. SARIMA 参数
     * =========================
     *
     * 小时级数据:
     *      1天 = 24点
     *
     * 所以季节周期:
     *      s = 24
     */
    int p = 3;
    int d = 1;
    int q = 1;

    int s = 24;
    int P = 3;
    int D = 1;
    int Q = 1;

    /*
     * =========================
     * 2. 读取数据
     * =========================
     */
    int Row_all, Col_all;
    const char *filename = "data/IPCLK.txt";

    double *X = load_matrix(filename, &Row_all, &Col_all);
    if (!X) {
        fprintf(stderr, "Failed to load matrix: %s\n", filename);
        return 1;
    }

    printf("rows = %d, cols = %d\n", Row_all, Col_all);

    /*
     * =========================
     * 3. 分配原始数据数组
     * =========================
     */
    double *error_data = (double*)malloc(sizeof(double) * Row_all);
    double *skew_data  = (double*)malloc(sizeof(double) * (Row_all - 1));
    double *t_data     = (double*)malloc(sizeof(double) * Row_all);

    if (!error_data || !skew_data || !t_data) {
        fprintf(stderr, "malloc failed for raw data arrays\n");
        free(X);
        return 1;
    }

    /*
     * =========================
     * 4. 提取鉴相值
     * =========================
     *
     * 如果你的 data/IPCLK.txt 是两列:
     *
     *      第1列: 序号
     *      第2列: 鉴相值
     *
     * 使用:
     *      phase_col = 1
     *
     * 如果你的 data/IPCLK.txt 是三列:
     *
     *      第1列: 序号
     *      第2列: 时间戳, 秒
     *      第3列: 鉴相值
     *
     * 改成:
     *      phase_col = 2
     *      time_col  = 1
     *
     * 当前默认按你原代码:
     *      第二列是鉴相值
     */
    int phase_col = 1;
    int use_time_col = 0;
    int time_col = 1;

    if (phase_col >= Col_all) {
        fprintf(stderr, "phase_col=%d out of range, Col_all=%d\n", phase_col, Col_all);
        free(X);
        free(error_data);
        free(skew_data);
        free(t_data);
        return 1;
    }

    for (i = 0; i < Row_all; ++i) {
        error_data[i] = X[i * Col_all + phase_col];

        if (use_time_col) {
            t_data[i] = X[i * Col_all + time_col];
        } else {
            /*
             * 默认用行号作为秒级时间。
             * 如果原始数据是1秒一个点，那么 i 就等价于秒。
             */
            t_data[i] = (double)i;
        }
    }

    /*
     * =========================
     * 5. 转换为相对鉴相值，单位 ns
     * =========================
     *
     * 原始鉴相值单位:
     *      1 count = 1 / 256 ns
     *
     * 转换:
     *      relative_phase_ns = (phase - phase0) / 256
     */
    double error0 = error_data[0];

    for (i = 0; i < Row_all; ++i) {
        error_data[i] = (error_data[i] - error0) / 256.0;
    }

    /*
     * =========================
     * 6. 一阶差分
     * =========================
     *
     * 这个数组当前主要保留，后续如果要做平均斜率模型可以用。
     */
    for (i = 0; i < Row_all - 1; ++i) {
        skew_data[i] = error_data[i + 1] - error_data[i];
    }

    /*
     * =========================
     * 7. 小时平均
     * =========================
     *
     * win_sz = 3600
     *
     * 如果原始数据是1秒一个点:
     *      3600点 = 1小时
     *
     * 得到:
     *      mean_data: 小时级相对鉴相值，单位 ns
     */
    int win_sz = 3600;
    int mean_data_n;

    double *mean_data = move_average(error_data, Row_all, win_sz, &mean_data_n);

    if (!mean_data) {
        fprintf(stderr, "move_average failed\n");
        free(X);
        free(error_data);
        free(skew_data);
        free(t_data);
        return 1;
    }

    printf("mean_data_n = %d\n", mean_data_n);

    /*
     * 保存小时平均数据
     */
    FILE *fp_mean = fopen("data/IPCLK_mean_data.txt", "w");
    if (!fp_mean) {
        perror("fopen data/IPCLK_mean_data.txt failed");
        free(X);
        free(error_data);
        free(skew_data);
        free(t_data);
        free(mean_data);
        return 1;
    }

    for (i = 0; i < mean_data_n; i++) {
        fprintf(fp_mean, "%.6f\n", mean_data[i]);
    }

    fclose(fp_mean);

    /*
     * =========================
     * 8. 构造小时级时间轴
     * =========================
     *
     * t_data_step:
     *      0, 3600, 7200, 10800, ...
     *
     * 单位:
     *      秒
     */
    double *t_data_step = (double*)malloc(sizeof(double) * mean_data_n);

    if (!t_data_step) {
        fprintf(stderr, "malloc failed for t_data_step\n");
        free(X);
        free(error_data);
        free(skew_data);
        free(t_data);
        free(mean_data);
        return 1;
    }

    for (i = 0; i < mean_data_n; ++i) {
        t_data_step[i] = t_data[i * win_sz];
    }

    /*
     * =========================
     * 9. 滚动窗口参数
     * =========================
     *
     * 当前数据已经是小时级:
     *      1天 = 24点
     *
     * 训练窗口:
     *      10天 = 240点
     *
     * 移动步长:
     *      0.25天 = 6小时 = 6点
     *
     * 起点:
     *      第0天 到 第12天
     */
    int N = mean_data_n;

    int points_per_day = 24;
    int train_days = 10;
    int N_train = train_days * points_per_day;

    double rolling_step_day = 0.25;
    int rolling_step_points = (int)(rolling_step_day * points_per_day + 0.5);

    int max_start_day = 12;
    int max_start_idx = max_start_day * points_per_day;

    printf("points_per_day = %d\n", points_per_day);
    printf("N_train = %d\n", N_train);
    printf("rolling_step_day = %.6f\n", rolling_step_day);
    printf("rolling_step_points = %d\n", rolling_step_points);
    printf("max_start_day = %d\n", max_start_day);
    printf("max_start_idx = %d\n", max_start_idx);

    if (rolling_step_points <= 0) {
        fprintf(stderr, "rolling_step_points <= 0\n");
        free(X);
        free(error_data);
        free(skew_data);
        free(t_data);
        free(mean_data);
        free(t_data_step);
        return 1;
    }

    if (N <= N_train) {
        fprintf(stderr, "Not enough data: N=%d, N_train=%d\n", N, N_train);
        free(X);
        free(error_data);
        free(skew_data);
        free(t_data);
        free(mean_data);
        free(t_data_step);
        return 1;
    }

    /*
     * =========================
     * 10. 打开输出文件
     * =========================
     */
    FILE *fp_summary = fopen("data/rolling_summary.csv", "w");
    if (!fp_summary) {
        perror("fopen data/rolling_summary.csv failed");
        free(X);
        free(error_data);
        free(skew_data);
        free(t_data);
        free(mean_data);
        free(t_data_step);
        return 1;
    }

    FILE *fp_detail = fopen("data/rolling_detail.csv", "w");
    if (!fp_detail) {
        perror("fopen data/rolling_detail.csv failed");
        fclose(fp_summary);
        free(X);
        free(error_data);
        free(skew_data);
        free(t_data);
        free(mean_data);
        free(t_data_step);
        return 1;
    }

    /*
     * rolling_summary.csv:
     *      每个窗口一行，保存误差统计
     */
    fprintf(
        fp_summary,
        "case_id,start_day,train_start_day,train_end_day,forecast_days,n_test,"
        "linear_mae,linear_rmse,linear_max_abs,linear_exceed_100,"
        "sarima_mae,sarima_rmse,sarima_max_abs,sarima_exceed_100,"
        "linear_sarima_mae,linear_sarima_rmse,linear_sarima_max_abs,linear_sarima_exceed_100\n"
    );

    /*
     * rolling_detail.csv:
     *      每个预测点一行，保存真实值、预测值、误差
     */
    fprintf(
        fp_detail,
        "case_id,start_day,forecast_point,forecast_day,"
        "actual,"
        "linear_pred,linear_error,"
        "sarima_pred,sarima_error,"
        "linear_sarima_pred,linear_sarima_error\n"
    );

    /*
     * =========================
     * 11. 滚动窗口主循环
     * =========================
     */
    int case_id = 0;

    for (
        int train_start_idx = 0;
        train_start_idx <= max_start_idx;
        train_start_idx += rolling_step_points
    ) {
        int train_end_idx = train_start_idx + N_train;

        /*
         * 训练窗口越界则跳过
         */
        if (train_end_idx >= N) {
            printf(
                "Skip start_idx=%d because train_end_idx=%d >= N=%d\n",
                train_start_idx,
                train_end_idx,
                N
            );
            continue;
        }

        int N_test = N - train_end_idx;
        int L = N_test;

        if (N_test <= 0) {
            printf("Skip start_idx=%d because N_test <= 0\n", train_start_idx);
            continue;
        }

        double start_day = train_start_idx / (double)points_per_day;
        double train_start_day = train_start_idx / (double)points_per_day;
        double train_end_day = train_end_idx / (double)points_per_day;
        double forecast_days = N_test / (double)points_per_day;

        printf("\n==============================\n");
        printf("case_id         = %d\n", case_id);
        printf("train_start_day = %.6f\n", train_start_day);
        printf("train_end_day   = %.6f\n", train_end_day);
        printf("forecast_days   = %.6f\n", forecast_days);
        printf("N_test          = %d\n", N_test);
        printf("==============================\n");

        /*
         * =========================
         * 11.1 分配当前窗口数组
         * =========================
         */
        double *sarima_pred = (double*)malloc(sizeof(double) * L);
        double *amse        = (double*)malloc(sizeof(double) * L);

        double *x_train = (double*)malloc(sizeof(double) * N_train);
        double *y_train = (double*)malloc(sizeof(double) * N_train);
        double *x_test  = (double*)malloc(sizeof(double) * N_test);
        double *y_test  = (double*)malloc(sizeof(double) * N_test);

        if (!sarima_pred || !amse || !x_train || !y_train || !x_test || !y_test) {
            fprintf(stderr, "malloc failed in case %d\n", case_id);

            free(sarima_pred);
            free(amse);
            free(x_train);
            free(y_train);
            free(x_test);
            free(y_test);

            fclose(fp_summary);
            fclose(fp_detail);

            free(X);
            free(error_data);
            free(skew_data);
            free(t_data);
            free(mean_data);
            free(t_data_step);

            return 1;
        }

        /*
         * 构造训练集
         */
        for (i = 0; i < N_train; ++i) {
            int idx = train_start_idx + i;
            x_train[i] = t_data_step[idx];
            y_train[i] = mean_data[idx];
        }

        /*
         * 构造测试集
         */
        for (i = 0; i < N_test; ++i) {
            int idx = train_end_idx + i;
            x_test[i] = t_data_step[idx];
            y_test[i] = mean_data[idx];
        }

        /*
         * =========================
         * 11.2 Linear 模型
         * =========================
         */
        double k, b;

        double *linear_pred_y_test = (double*)malloc(sizeof(double) * N_test);

        if (!linear_pred_y_test) {
            fprintf(stderr, "malloc linear_pred_y_test failed in case %d\n", case_id);

            free(sarima_pred);
            free(amse);
            free(x_train);
            free(y_train);
            free(x_test);
            free(y_test);

            fclose(fp_summary);
            fclose(fp_detail);

            free(X);
            free(error_data);
            free(skew_data);
            free(t_data);
            free(mean_data);
            free(t_data_step);

            return 1;
        }

        linear_fit(x_train, y_train, N_train, &k, &b);
        linear_predict(k, b, x_test, linear_pred_y_test, N_test);

        /*
         * =========================
         * 11.3 SARIMA 模型
         * =========================
         */
        sarima_object obj = sarima_init(p, d, q, s, P, D, Q, N_train);

        if (!obj) {
            fprintf(stderr, "sarima_init failed in case %d\n", case_id);

            free(sarima_pred);
            free(amse);
            free(x_train);
            free(y_train);
            free(x_test);
            free(y_test);
            free(linear_pred_y_test);

            fclose(fp_summary);
            fclose(fp_detail);

            free(X);
            free(error_data);
            free(skew_data);
            free(t_data);
            free(mean_data);
            free(t_data_step);

            return 1;
        }

        /*
         * Method 0 = MLE
         */
        sarima_setMethod(obj, 0);

        sarima_exec(obj, y_train);

        /*
         * 预测未来 L 个点
         */
        sarima_predict(obj, y_train, L, sarima_pred, amse);

        /*
         * =========================
         * 11.4 Linear_SARIMA 融合模型
         * =========================
         *
         * 思路:
         *      1. 对 SARIMA 预测结果拟合一条旧线性趋势
         *      2. 用 sarima_pred - old_linear_trend 提取 SARIMA 的残差形状
         *      3. 用训练窗口后半段拟合新的线性斜率
         *      4. 新线性趋势 + SARIMA残差形状 + 起点偏移
         */
        double k_old, b_old;

        double *residual = (double*)malloc(sizeof(double) * N_test);
        double *linear_sarima_pred = (double*)malloc(sizeof(double) * N_test);

        if (!residual || !linear_sarima_pred) {
            fprintf(stderr, "malloc residual or linear_sarima_pred failed in case %d\n", case_id);

            sarima_free(obj);

            free(sarima_pred);
            free(amse);
            free(x_train);
            free(y_train);
            free(x_test);
            free(y_test);
            free(linear_pred_y_test);
            free(residual);
            free(linear_sarima_pred);

            fclose(fp_summary);
            fclose(fp_detail);

            free(X);
            free(error_data);
            free(skew_data);
            free(t_data);
            free(mean_data);
            free(t_data_step);

            return 1;
        }

        /*
         * 1. 对 SARIMA 预测结果拟合旧趋势
         */
        linear_fit(x_test, sarima_pred, N_test, &k_old, &b_old);

        /*
         * 2. 提取 SARIMA 预测中的残差形状
         */
        for (i = 0; i < N_test; i++) {
            double y_trend_old = k_old * x_test[i] + b_old;
            residual[i] = sarima_pred[i] - y_trend_old;
        }

        /*
         * 3. 用训练窗口后半段拟合新的线性趋势
         *
         * unstable_day = 5 表示:
         *      对每个10天训练窗口，跳过前5天，用后5天估计新斜率。
         *
         * 注意:
         *      这里使用的是窗口内部索引，不再使用绝对天数。
         *      这样滚动窗口从第几天开始都不会出错。
         */
        double unstable_day = 5.0;
        int fit_start_idx = (int)(unstable_day * points_per_day + 0.5);

        if (fit_start_idx < 0) {
            fit_start_idx = 0;
        }

        if (fit_start_idx >= N_train - 2) {
            fit_start_idx = N_train / 2;
        }

        int fit_n = N_train - fit_start_idx;

        if (fit_n < 2) {
            fprintf(stderr, "fit_n < 2 in case %d\n", case_id);

            sarima_free(obj);

            free(sarima_pred);
            free(amse);
            free(x_train);
            free(y_train);
            free(x_test);
            free(y_test);
            free(linear_pred_y_test);
            free(residual);
            free(linear_sarima_pred);

            fclose(fp_summary);
            fclose(fp_detail);

            free(X);
            free(error_data);
            free(skew_data);
            free(t_data);
            free(mean_data);
            free(t_data_step);

            return 1;
        }

        double k_new, b_tmp;

        linear_fit(
            x_train + fit_start_idx,
            y_train + fit_start_idx,
            fit_n,
            &k_new,
            &b_tmp
        );

        /*
         * 4. 起点偏移
         *
         * 使 Linear_SARIMA 的第一个预测点和 SARIMA 的第一个预测点对齐。
         */
        double base = sarima_pred[0] - k_new * x_test[0];

        for (i = 0; i < N_test; i++) {
            double y_trend_new = k_new * x_test[i];
            linear_sarima_pred[i] = y_trend_new + residual[i] + base;
        }

        /*
         * =========================
         * 11.5 计算误差指标
         * =========================
         */
        double linear_mae, linear_rmse, linear_max_abs;
        int linear_exceed_100;

        double sarima_mae, sarima_rmse, sarima_max_abs;
        int sarima_exceed_100;

        double linear_sarima_mae, linear_sarima_rmse, linear_sarima_max_abs;
        int linear_sarima_exceed_100;

        calc_error_metrics(
            y_test,
            linear_pred_y_test,
            N_test,
            &linear_mae,
            &linear_rmse,
            &linear_max_abs,
            &linear_exceed_100
        );

        calc_error_metrics(
            y_test,
            sarima_pred,
            N_test,
            &sarima_mae,
            &sarima_rmse,
            &sarima_max_abs,
            &sarima_exceed_100
        );

        calc_error_metrics(
            y_test,
            linear_sarima_pred,
            N_test,
            &linear_sarima_mae,
            &linear_sarima_rmse,
            &linear_sarima_max_abs,
            &linear_sarima_exceed_100
        );

        /*
         * =========================
         * 11.6 保存汇总结果
         * =========================
         */
        fprintf(
            fp_summary,
            "%d,%.6f,%.6f,%.6f,%.6f,%d,"
            "%.6f,%.6f,%.6f,%d,"
            "%.6f,%.6f,%.6f,%d,"
            "%.6f,%.6f,%.6f,%d\n",

            case_id,
            start_day,
            train_start_day,
            train_end_day,
            forecast_days,
            N_test,

            linear_mae,
            linear_rmse,
            linear_max_abs,
            linear_exceed_100,

            sarima_mae,
            sarima_rmse,
            sarima_max_abs,
            sarima_exceed_100,

            linear_sarima_mae,
            linear_sarima_rmse,
            linear_sarima_max_abs,
            linear_sarima_exceed_100
        );

        /*
         * =========================
         * 11.7 保存每个预测点的详细结果
         * =========================
         */
        for (i = 0; i < N_test; i++) {
            double forecast_day = (train_end_idx + i) / (double)points_per_day;

            double linear_error = y_test[i] - linear_pred_y_test[i];
            double sarima_error = y_test[i] - sarima_pred[i];
            double linear_sarima_error = y_test[i] - linear_sarima_pred[i];

            fprintf(
                fp_detail,
                "%d,%.6f,%d,%.6f,"
                "%.6f,"
                "%.6f,%.6f,"
                "%.6f,%.6f,"
                "%.6f,%.6f\n",

                case_id,
                start_day,
                i,
                forecast_day,

                y_test[i],

                linear_pred_y_test[i],
                linear_error,

                sarima_pred[i],
                sarima_error,

                linear_sarima_pred[i],
                linear_sarima_error
            );
        }

        /*
         * =========================
         * 11.8 打印当前窗口结果
         * =========================
         */
        printf(
            "Linear        MAE = %.6f, RMSE = %.6f, MaxAbs = %.6f, exceed100 = %d\n",
            linear_mae,
            linear_rmse,
            linear_max_abs,
            linear_exceed_100
        );

        printf(
            "SARIMA        MAE = %.6f, RMSE = %.6f, MaxAbs = %.6f, exceed100 = %d\n",
            sarima_mae,
            sarima_rmse,
            sarima_max_abs,
            sarima_exceed_100
        );

        printf(
            "Linear_SARIMA MAE = %.6f, RMSE = %.6f, MaxAbs = %.6f, exceed100 = %d\n",
            linear_sarima_mae,
            linear_sarima_rmse,
            linear_sarima_max_abs,
            linear_sarima_exceed_100
        );

        /*
         * =========================
         * 11.9 释放当前窗口内存
         * =========================
         */
        sarima_free(obj);

        free(sarima_pred);
        free(amse);

        free(x_train);
        free(y_train);
        free(x_test);
        free(y_test);

        free(linear_pred_y_test);
        free(residual);
        free(linear_sarima_pred);

        case_id++;
    }

    /*
     * =========================
     * 12. 关闭文件
     * =========================
     */
    fclose(fp_summary);
    fclose(fp_detail);

    /*
     * =========================
     * 13. 释放全局内存
     * =========================
     */
    free(X);
    free(error_data);
    free(skew_data);
    free(t_data);
    free(mean_data);
    free(t_data_step);

    printf("\nRolling window test finished.\n");
    printf("Summary saved to data/rolling_summary.csv\n");
    printf("Detail saved to data/rolling_detail.csv\n");

    return 0;
}