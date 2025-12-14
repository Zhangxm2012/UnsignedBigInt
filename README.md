# UnsignedBigInt
无符号大整数库
（还不完善，加减乘除模、左右移、gcd、lcm、阶乘已实现）
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
