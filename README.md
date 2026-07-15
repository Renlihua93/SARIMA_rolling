# CTSA SARIMA (simplified)

精简版 seasonal SARIMA，从 `ctsa-master` 重组。对外 API：`include/sarima.h`。

## 源文件（`src/`，6 个 `.c`）

| 文件 | 说明 |
|------|------|
| `sarima.c` | 对外 API |
| `fit.c` | 季节 ARIMA 拟合 / Kalman 预测 |
| `ts.c` | 差分、卷积、AR/MA 根检验 |
| `optim.c` | BFGS、`fminunc`、`hessian_fd` |
| `linalg.c` | 矩阵乘、转置、LU、求逆 |
| `cstats.c` | `mean`、`polyroot` |

从上游重新裁剪（需本机存在 `../../ctsa-master`）：

```powershell
powershell -File scripts/extract_sarima.ps1
```

## 构建

CLion 打开本目录 → Reload CMake → 构建/运行 **`sarimatest`**。

```text
cmake --build cmake-build-debug --target sarimatest
```

勿对 `sarimatest.c` 单文件编译（需链接 `sarimalib`）。
# SARIMA_simple
