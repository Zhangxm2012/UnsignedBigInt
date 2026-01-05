# 使用指南

项目传送门：[github](https://github.com/Zhangxm2012/UnsignedBigInt)，[luogu](https://www.luogu.com.cn/article/hjafe8xc)，欢迎支持。

## 异常类

在 `ERROR` 中，放了几个异常类，分别为：

- 基底：`Exception`；
- 除零：`Div_by_zero`；
- 超出内存限制：`MLE`；
- 负数：`Negative`；
- 非数字：`Number`；
- 越界：`Out_of_range`。

## 变换

`Transform` 中分别有：有 AVX2 优化和朴素实现。

### AVX2 优化

实现了一个短小但功能齐全的优化复数类 `Complex`，此部分不过多赘述。

重头戏是 `FFT` 类：

- `init(int Len)`：预处理出 $2^{\lfloor\log_2Len\rfloor}$ 内的单位根（旋转因子）；
- `dif(vector<Complex>&a)`：基 $2$ 频域抽取；
- `dit(vector<Complex>&a)`：基 $2$ 时域抽取；
- `mul(vector<Complex>&a,vector<Complex>&b)`：频域相乘。

### 朴素

实现了一个短小的复数类 `fComplex`。

`FFT` 类：

- `init(int len)`：分配 `fft_a` 的内存；
- `Init(int k)`：预处理出 $2^k$ 内的单位根（旋转因子）；
- `fft(int len,int flag)`：若 `flag` 为 $1$，做 DFT；否则，做 IDFT；

## 输入输出

`IO` 中的快读快输没什么好说的，主要是为了卡常；

## 大整数库

在 `UnsignedBigInt` 中，封装了赋值、加减乘除模（支持高精乘模除低精、左右移），但位运算还未实现（$O(n^2)$ 太抽象了）。\
值得注意的是，这里的长度是在 `BASE` 进制下。

### 宏定义

- `lf`：`double`；
- `ull`：`unsigned long long`；
- `u32`：`unsigned`；
- `int128`：`__int128_t`；
- `__AVX2__`：调试用，记得配上 `#pragma GCC target("fma")`；
- `LENGTH`：最大长度，可通过在引用文件前定义 `SIZE` 修改。

### 私有变量

- `LEN=8`，压位长度；
- `BASE=100000000`，进制；
- `FFT_BASE=10000`，FFT 中的进制；
- `T=255`，算法调节阈值，长度大于 `T`，使用 FFT 优化 / 牛顿迭代 / 快速幂左右移；
- `INV=64`，牛顿迭代下限，若长度小于 `INV`，使用暴力除法；
- `DEFAULT=1000`，默认分配长度；
- `MAX=ULLONG_MAX/(BASE+1)`，高精乘低精阈值，小于等于其直接乘；否则转成高精乘高精；
- `LOG=std::__lg(MAX)-1;`，分块式左移的长度（？）；
- `std::unique_ptr<ull[]>num;`，数码；
- `int len,Max;`，分别为该数的长度，分配的最大长度；
- `static constexpr ull p10[]={1,10,100,1000,10000,100000,1000000,10000000,100000000};`，返回 $10$ 的幂；

### 私有函数

- `Pow_10(int p)`：返回 $10^p,0\le p\le 8$；
- `is_zero()`：返回该数是否为 $0$；
- `Expand(int nMax)`：分配至少 `nMax` 的内存（若 `Max` 大于 `nMax`，啥也不干），**不清零**；
- `Mul(const UnsignedBigInt&b)`：暴力高精乘；
- `Get(const UnsignedBigInt&a,int pos)`：除法辅助函数，用于估商；
- `Simple_Mod(const UnsignedBigInt&b)`：暴力除法，返回值为 `std::pair<UnsignedBigInt,UnsignedBigInt>`；
- `Left(int cnt)`：返回该数乘 `BASE` 的 `cnt` 次方；
- `Right(int cnt)`：返回该数除以 `BASE` 的 `cnt` 次方；
- `Inv(int n)`：返回该数的倒数，保留 `n` 位精度；
- `Mod(const UnsignedBigInt&b)`：优化除法，返回值与 `Simple_Mod` 相同，[record](https://judge.yosupo.jp/submission/343251)；
- `Get_sqrt()`：开平方初值估计。

### 公开部分

- `UnsignedBigInt()`：构造函数（置零）；
- `~UnsignedBigInt()`：解析函数；
- `UnsignedBigInt(const UnsignedBigInt&b)`：构造函数（赋值）；
- `UnsignedBigInt(UnsignedBigInt&&b)noexcept`：构造函数（移动）；
- `UnsignedBigInt(const ull&b)`：构造函数（低精）；
- `init(const u32&size=DEFAULT)`：初始化函数；
- `UnsignedBigInt(const std::string&s)`：构造函数（string）；
- `UnsignedBigInt(const char*s)`：构造函数（char\*）;
- `UnsignedBigInt(std::string_view s)`：构造函数（C++17 以上）；
- `fread()`：使用 `IO` 快速读入；
- `fwrite()`：使用 `IO` 快速输出；
- `friend std::istream& operator>>(std::istream&in,UnsignedBigInt&x)`：使用 `>>` 输入；
- `friend std::ostream& operator<<(std::ostream&out,const UnsignedBigInt&x)`：使用 `<<` 输出；
- `Two()`：返回该数 $2$ 因子的个数；
- `Cmp(const UnsignedBigInt&b)`：类似于三路比较符；
- `at(const int&b)`：获取 `num[b]`，**进行越界检查**；
- `pow(const ull&b)`：返回该数的 `b` 次方；
- `root(int m)`：返回该数的 `b` 次根；
- `sqrt()`：返回该数的平方根（优于 `root(2)`）；
- `True()`：返回该数的布尔类型（因实现了高精乘低精，`bool operator()` 会有问题）；
- `operator std::string()`：返回该数的字符串类型；
- `square()`：自乘，快于 `*this*=*this`；
- `Square()`：返回该数的平方（不修改）；
- `operator""_UI(const char*literal,size_t len)`：字面量；

## 运算

`Operation` 中：

- `Pow(const UnsignedBigInt&a,int p)`：返回 $a^p$ 的准确值；
- `Pow(const UnsignedBigInt&a,int p,const UnsignedBigInt&Mod)`：返回 $a^p \bmod Mod$ 的值；
- `Fact(int st,int n)`：返回从 $st$ 乘到 $n$；
- `Fact(int n)`：返回 $n!$；
- `Gcd(const UnsignedBigInt&a,const UnsignedBigInt&b)`：返回 $\gcd(a,b)$；
- `Lcm(const UnsignedBigInt&a,const UnsignedBigInt&b)`：返回 $\operatorname{lcm}(a,b)$；
- `Root(const UnsignedBigInt&a,ull p)`：返回 $\sqrt[p]{a}$；
- `Random(int len)`：随机返回一个长度（十进制下，即数码个数）为 `len` 的大整数；
- `Sqrt(const UnsignedBigInt&a)`：返回 $\sqrt{a}$；

# 关于效率

加法：未测评（应该挂不了）；  
乘法：[record](https://judge.yosupo.jp/submission/343248)；  
除法：[record](https://judge.yosupo.jp/submission/343251)。
