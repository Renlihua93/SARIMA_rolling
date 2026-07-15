#ifndef MATRIX_IO_H
#define MATRIX_IO_H

#include <stdio.h>
#include <stdlib.h>

/*
 读取二维数组 txt 文件
 返回连续内存：row-major 顺序
 data[i * cols + j]
*/
double *load_matrix(
    const char *filename,
    int *rows,
    int *cols
);

/*
 释放矩阵内存
*/
void free_matrix(double *matrix);

/*
 打印矩阵（调试用）
*/
void print_matrix(const double *data, int rows, int cols);

double *move_average(const double *data, int len, int win_sz, int *out_n);

#endif