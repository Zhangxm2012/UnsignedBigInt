# UnsignedBigInt/SignedBigInt
无符号大整数库 / 有符号大整数库（作者偷懒放一起了）
（还不完善，位运算还未实现）
## upd on 2025/11/23
添加了快读快输出
## upd on 2025/11/25
改进了快速幂代码，使其能通过 [这道题](https://www.luogu.com.cn/problem/U393978)
## upd on 2025/11/26
添加了阶乘，[测试记录](https://vjudge.net/solution/65906519)，还有 gcd，在 [这道题](https://www.luogu.com.cn/problem/P2152) 中，最慢跑了 $96$ms
## upd on 2025/12/8
添加了朴素除法
## upd on 2025/12/9
补了高精除低精的锅，添加了真正的大数除法，并在 [这道题](https://loj.ac/p/164) 中最慢跑了 $50$ms
## upd on 2025/12/14
更新左右移，和通过在引用文件前定义宏 SIZE，来修改长度
```cpp
#define SIZE 114514
#include "UnsignedBigInt.cpp"
//your code
```
## upd on 2025/12/15
实现基本完成，下一步是换成动态数组和支持负数
## upd on 2025/12/23
动态已实现，请在 New 中查看最新版本（Prime 为素数判断实现，目前为空）
## upd on 2026/01/01
实现了 FFT 的优化，欸，挂了
## upd on 2026/01/05
FFT 活了
