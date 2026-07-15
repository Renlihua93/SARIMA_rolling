#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* 统计行数 */
static int count_rows(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("fopen");
        return 0;
    }

    int rows = 0;
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        rows++;
    }

    fclose(fp);
    return rows;
}

/* 统计列数（根据第一行） */
static int count_cols(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("fopen");
        return 0;
    }

    char line[1024];
    fgets(line, sizeof(line), fp);

    int cols = 0;
    char *token = strtok(line, " \t\n");
    while (token) {
        cols++;
        token = strtok(NULL, " \t\n");
    }

    fclose(fp);
    return cols;
}

/*
 读取二维数组 txt 文件
 返回连续内存：data[i * cols + j]
*/
double *load_matrix(const char *filename, int *rows, int *cols) {
    *rows = count_rows(filename);
    *cols = count_cols(filename);

    if (*rows == 0 || *cols == 0) {
        return NULL;
    }

    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("fopen");
        return NULL;
    }

    double *data = (double *)malloc((*rows) * (*cols) * sizeof(double));
    if (!data) {
        fclose(fp);
        return NULL;
    }

    for (int i = 0; i < *rows; i++) {
        for (int j = 0; j < *cols; j++) {
            if (fscanf(fp, "%lf", &data[i * (*cols) + j]) != 1) {
                fprintf(stderr, "Read error at row %d col %d\n", i, j);
                free(data);
                fclose(fp);
                return NULL;
            }
        }
    }

    fclose(fp);
    return data;
}

/* 打印矩阵（调试用） */
void print_matrix(const double *data, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%10.4f ", data[i * cols + j]);
        }
        printf("\n");
    }
}


double *move_average(const double *data, int len, int win_sz, int *out_n) {
    if (!data || win_sz <= 0 || len <= 0) {
        *out_n = 0;
        return NULL;
    }

    int n = (int)ceil((double)len / win_sz);
    double *mean_data = (double *)malloc(n * sizeof(double));
    if (!mean_data) {
        *out_n = 0;
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        int start = i * win_sz;
        int end = start + win_sz;
        if (end > len) {
            end = len;
        }

        double sum = 0.0;
        int count = 0;

        for (int k = start; k < end; k++) {
            sum += data[k];
            count++;
        }

        mean_data[i] = sum / count;
    }

    *out_n = n;
    return mean_data;
}
