# CLion：ctsa-sarima-simple

1. **File → Open** 选择 `C:\project\ctsa-master\ctsa-sarima-simple`
2. **Reload CMake Project**
3. 运行配置选 **sarimatest**（CMake Application），**不要** Run `sarimatest.c`
4. **Build Project** (`Ctrl+F9`) → **Run**

库目标为 `sarimalib`（不是 `ctsalib`）。头文件：`include/sarima.h`。

若出现 `undefined reference to sarima_*`，说明用了单文件编译，请改选 CMake 目标。
