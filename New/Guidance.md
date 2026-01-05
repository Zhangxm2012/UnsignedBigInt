# 使用指南

## 异常类

在 `namespace ERROR` 中，放了几个异常类，分别为：

+ 基底：`Exception`；
+ 除零：`Div_by_zero`；
+ 超出内存限制：`MLE`；
+ 负数：`Negative`；
+ 非数字：`Number`；
+ 越界：`Out_of_range`。

## 变换

分为有 AVX2 优化和朴素实现。

### AVX2

实现了一个短小但功能齐全的优化复数类 `Complex`，此部分不过多赘述。

重头戏是 `FFT` 类：
+ `init(int Len)`：预处理出 $$