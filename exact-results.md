# 局部极小点计数的重采样协方差与共同邻居分布

日期：2026-10-10。状态：本次构造的候选命题，已给出自足证明和独立精确计算核验；历史原创性尚未认证。

## 结论与原创性边界

本次提出并证明的核心陈述是：**对于非空、有限、简单、无向、无三角形的正则图，且度数至少为 2，仅从局部极小点数量的整条重采样协方差曲线，就能唯一恢复顶点数、度数，以及不相邻顶点对的共同邻居数量分布。**

这里恢复的是一个明确的结构统计量，不是整张图。它不是把局部极小点模型改名：这个模型至少在 Baldi–Rinott–Stein 的 1988/1989 年工作中已出现。待确认的新颖点是“完整重采样曲线与共同邻居分布的显式可逆对应”，连同逐系数的精确公式。

本次 GitHub、网页及 arXiv 定向检索没有发现同一陈述；这不等价于证明历史上无人提出过，也不构成发表优先权认证。文末列出最接近文献，另附检索记录。以下论证不以文献未命中作为数学前提。

## 1. 定义与候选猜想

设 $`G=(V,E)`$，$`|V|=n\gt 0`$，每个顶点度数均为 $`d\ge2`$，且图中没有三角形。不要求连通。$`N(v)`$ 表示开邻域，$`N[v]=N(v)\cup\{v\}`$。

给每个顶点独立抽取 $`U_v\sim\mathrm{Unif}(0,1)`$，定义

```math
I_v(U)=\mathbf1\{U_v\lt U_w\text{ 对所有 }w\in N(v)\},
\qquad X_G(U)=\sum_{v\in V}I_v(U).
```

即 $`X_G`$ 是局部极小点数量。连续分布保证平局概率为零。这是一次局部比较，不是迭代贪心选取独立集。

对 $`0\le\rho\le1`$，另取彼此独立的
$`B_v\sim\mathrm{Bernoulli}(\rho)`$ 与 $`Z_v\sim\mathrm{Unif}(0,1)`$，并与所有 $`U_v`$ 独立。令

```math
U_v^{(\rho)}=\begin{cases}U_v,&B_v=1,\\Z_v,&B_v=0.\end{cases}
```

定义协方差曲线

```math
C_G(\rho)=\mathrm{Cov}\bigl(X_G(U),X_G(U^{(\rho)})\bigr).
```

因此 $`\rho=0`$ 表示完全重采样，$`\rho=1`$ 表示完全保留；$`C_G(0)=0`$，$`C_G(1)=\mathrm{Var}(X_G)`$。

对不相邻无序顶点对 $`\{u,v\}`$，令

```math
c_{uv}=|N(u)\cap N(v)|,\qquad
m_j=\#\{\{u,v\}:u\ne v,\ uv\notin E,\ c_{uv}=j\}.
```

**候选猜想。** $`C_G`$ 唯一决定 $`(n,d,m_0,\ldots,m_d)`$，反之这组数据也唯一决定 $`C_G`$。

下面直接证明它，因此在数学逻辑层面可将它称为命题；“候选”保留的是原创性状态，并非用计算替代证明。

## 2. 精确公式

约定 $`\binom ak=0`$ 当 $`k\gt a`$，空和为零，$`[\rho^k]P`$ 表示多项式的系数。定义

```math
S_k(G)=\sum_{\substack{\{u,v\}\subset V\\uv\notin E}}\binom{c_{uv}}k
=\sum_{j=k}^d\binom jk m_j.
```

设

```math
B_{n,d}(\rho)=
\frac{nd(d-1)}{2(d+1)^2(2d+1)}\rho^2
+\frac{nd}{(d+1)^2}
\sum_{k=3}^{d+1}\frac{\binom{d+1}k}{\binom{2d+1}k}\rho^k.
\tag{1}
```

**精确恒等式：**

```math
\boxed{
C_G(\rho)=B_{n,d}(\rho)+
\frac2{(d+1)^2}
\sum_{k=2}^{d}\frac{S_k(G)}{\binom{2d+1}k}\rho^k.
}
\tag{2}
```

以下给出全部推导。

## 3. 一个有限二项式恒等式

若 $`K\sim\mathrm{Bin}(r,\rho)`$，整数 $`a\gt r\ge0`$，则

```math
\mathbb E\frac K{a-K}
=\sum_{k=1}^{r}\frac{\binom rk}{\binom{a-1}k}\rho^k.
\tag{3}
```

证明：由 $`K/(a-K)=a/(a-K)-1`$，

```math
\begin{aligned}
\mathbb E\frac K{a-K}
&=a\int_0^1 x^{a-1}(1-\rho+\rho/x)^r\,dx-1\\
&=a\sum_{k=1}^r\binom rk\rho^k
\int_0^1x^{a-1-k}(1-x)^k\,dx.
\end{aligned}
```

所有和都是有限和；$`a\gt r`$ 保证积分存在。反复分部积分给出
$`\int_0^1x^{a-1-k}(1-x)^k dx=(a-1-k)!k!/a!`$，代入即得 (3)。证毕。

## 4. 三类顶点对的协方差

写 $`q=d+1`$，于是 $`\mathbb E I_v=1/q`$。把
$`\mathrm{Cov}(I_u(U),I_v(U^{(\rho)}))`$ 分成三类。

### 4.1 不相邻的不同顶点

固定不相邻的 $`u,v`$，共同邻居数为 $`c`$。两次事件共享的标签变量，只有那些被保留的共同邻居；其数量 $`K\sim\mathrm{Bin}(c,\rho)`$。

条件于 $`K=r`$，设两个中心标签为 $`x,y`$。在 $`x\lt y`$ 区域，两个局部极小事件同时发生的概率密度因子为
$`(1-x)^{d-r}(1-y)^d`$。由交换对称性，联合概率为

```math
\begin{aligned}
P_r
&=2\int_0^1\int_0^y(1-x)^{d-r}(1-y)^d\,dx\,dy\\
&=\frac2{q-r}\left(\frac1q-\frac1{2q-r}\right)
=\frac2{q(2q-r)}.
\end{aligned}
\tag{4}
```

所以 $`P_r-q^{-2}=q^{-2}r/(2q-r)`$。利用 (3)，得到

```math
N_c(\rho):=\mathrm{Cov}(I_u,I_v^{(\rho)})
=\frac1{q^2}\sum_{k=1}^{c}
\frac{\binom ck}{\binom{2d+1}k}\rho^k.
\tag{5}
```

### 4.2 同一个顶点

只需考察 $`N[v]`$ 中的 $`q`$ 个变量。保留数量 $`K\sim\mathrm{Bin}(q,\rho)`$。给定 $`K=r`$，保留集合在所有 $`r`$ 元子集上均匀。

中心标签以概率 $`r/q`$ 被保留。若保留中心，两次事件一共要求 $`2d-r+1`$ 个不同邻居标签大于同一个中心，因此联合概率为 $`1/(2q-r)`$。若不保留中心，联合概率由 (4) 给出 $`2/[q(2q-r)]`$。在相应情形概率为零时忽略该项。

故条件联合概率减去 $`q^{-2}`$ 等于

```math
\frac rq\frac1{2q-r}
+\frac{q-r}q\frac2{q(2q-r)}-\frac1{q^2}
=\frac d{q^2}\frac r{2q-r}.
```

再用 (3)，

```math
D_d(\rho):=\mathrm{Cov}(I_v,I_v^{(\rho)})
=\frac d{q^2}\sum_{k=1}^{d+1}
\frac{\binom{d+1}k}{\binom{2d+1}k}\rho^k.
\tag{6}
```

### 4.3 相邻顶点

若 $`uv\in E`$，无三角形条件保证 $`N[u]\cap N[v]=\{u,v\}`$。因此两次局部事件可能共享的变量仅为这两个端点标签。

两者都不保留时，联合概率为 $`q^{-2}`$。两者都保留时，两个事件分别要求 $`U_u\lt U_v`$ 与 $`U_v\lt U_u`$，联合概率为零。恰好保留一个端点时，联合概率为

```math
J=\int_0^1(1-x)^d\int_0^x(1-y)^{d-1}\,dy\,dx
=\frac1d\left(\frac1{d+1}-\frac1{2d+1}\right)
=\frac1{q(2d+1)}.
```

于是

```math
\begin{aligned}
A_d(\rho)&:=(1-\rho)^2q^{-2}+2\rho(1-\rho)J-q^{-2}\\
&=-\frac{2d\rho+\rho^2}{q^2(2d+1)}.
\end{aligned}
\tag{7}
```

这也明确指出无三角形假设用在何处，不能直接把 (7) 套到一般图上。

## 5. 汇总与一次项抵消

方差/协方差展开使用的是有序顶点对。对角项有 $`n`$ 个，有序相邻对有 $`nd`$ 个，每个非邻无序对出现两次。因此

```math
C_G(\rho)=nD_d(\rho)+ndA_d(\rho)
+2\sum_{\{u,v\}:uv\notin E}N_{c_{uv}}(\rho).
\tag{8}
```

按中间顶点计数长度为 2 的路径。因为没有三角形，两端必不相邻，故

```math
S_1=\sum_{uv\notin E}c_{uv}=n\binom d2.
\tag{9}
```

(8) 的 $`\rho`$ 系数乘以 $`q^2(2d+1)`$ 后为

```math
ndq-2nd^2+nd(d-1)=0.
```

二次项中与图形状无关的部分为

```math
\frac{nd}{q^2}\frac{\binom q2}{\binom{2d+1}2}
-\frac{nd}{q^2(2d+1)}
=\frac{nd(d-1)}{2q^2(2d+1)}.
```

其余系数直接由 (5)、(6) 得到，恰为 (1)、(2)。这证明了对所有允许规模和所有 $`\rho\in[0,1]`$ 的恒等式；有限计算没有参与本证明。

## 6. 从曲线恢复结构

**第一步：恢复 $`d`$ 和 $`n`$。** (2) 中结构相关的和最高次数为 $`d`$，而

```math
[\rho^{d+1}]C_G=
\frac{nd}{(d+1)^2\binom{2d+1}{d+1}}\gt 0.
```

所以

```math
d=\deg C_G-1,
\qquad
n=\frac{(d+1)^2}{d}\binom{2d+1}{d+1}
[\rho^{d+1}]C_G.
\tag{10}
```

这解释了为什么无需预先告知 $`n,d`$。同时也解释为何主命题排除 $`d=1`$：匹配上的局部极小点数恒定，$`C_G\equiv0`$，无法从它恢复顶点数。

**第二步：恢复二项矩。** 对 $`2\le k\le d`$，

```math
S_k=\frac{(d+1)^2}2\binom{2d+1}k
[\rho^k]\bigl(C_G-B_{n,d}\bigr).
\tag{11}
```

另外

```math
S_0=\binom n2-\frac{nd}2,\qquad
S_1=n\binom d2.
\tag{12}
```

**第三步：恢复共同邻居分布。** 由
$`S_k=\sum_{j\ge k}\binom jk m_j`$，可按 $`j=d,d-1,\ldots,0`$ 递推

```math
m_j=S_j-\sum_{\ell=j+1}^d\binom\ell j m_\ell.
\tag{13}
```

这是一组对角系数为 1 的三角方程，因此解唯一。反方向由 (2) 立即成立。候选猜想得证。$`\square`$

## 7. 两个直接推论

### 四环由二次系数精确给出

每个四环有两对对角顶点，故 $`S_2=2\,#C_4(G)`$。由 (2)，

```math
\boxed{
[\rho^2]C_G=
\frac{nd(d-1)}{2(d+1)^2(2d+1)}
+\frac{4\,#C_4(G)}{d(d+1)^2(2d+1)}.
}
\tag{14}
```

对 $`k\ge3`$，$`S_k`$ 就是 $`K_{2,k}`$ 子图的个数：从不相邻的一对顶点及其 $`k`$ 个共同邻居得到一个这样的子图；两侧大小不同，无四环的双重表示问题。

### 对所有保留概率同时成立的极值界

令

```math
U_{n,d}(\rho)=B_{n,d}(\rho)+
\frac{n(d-1)}{(d+1)^2}
\sum_{k=2}^d\frac{\binom dk}{\binom{2d+1}k}\rho^k.
```

则

```math
B_{n,d}(\rho)\le C_G(\rho)\le U_{n,d}(\rho).
\tag{15}
```

证明：$`S_k\ge0`$。对 $`0\le c\le d`$，
$`\binom ck\le(c/d)\binom dk`$：在 $`c\ge1`$ 时用
$`\binom ck/c=\binom{c-1}{k-1}/k`$ 即可。由 (9)，
$`S_k\le n(d-1)\binom dk/2`$。代回 (2) 得证。

对任意固定 $`\rho\gt 0`$，下等号当且仅当图无四环。上等号当且仅当图是若干 $`K_{d,d}`$ 的不交并。后者的必要性如下：二次系数取等要求所有非邻对的共同邻居数只取 $`0,d`$。固定顶点 $`v`$ 及 $`a\in N(v)`$；对每个 $`x\in N(a)\setminus\{v\}`$，$`x`$ 与 $`v`$ 不相邻但共有邻居 $`a`$，故 $`N(x)=N(v)`$。集合 $`N(a)`$ 与 $`N(v)`$ 因无三角形而不交，各有 $`d`$ 个顶点，彼此完全连接，度数已用满，构成一个 $`K_{d,d}`$ 连通分支。对各分支重复即可。充分性直接代入。

这些是统一上、下界；对指定的 $`n,d`$，取等图可能不存在，例如上界取等要求 $`2d\mid n`$。

## 8. 例子与可恢复信息的限度

$`C_8`$ 与两个不相交四环 $`2C_4`$ 都有 $`n=8,d=2`$，但

```math
C_{C_8}(\rho)=\frac8{45}\rho^2+\frac8{45}\rho^3,
\qquad
C_{2C_4}(\rho)=\frac4{15}\rho^2+\frac8{45}\rho^3.
```

两条曲线的二次系数之差 $`4/45`$，精确对应后者多出的两个四环。

另一方面，$`C_{10}`$ 与 $`2C_5`$ 的数据均为
$`n=10,d=2,m_0=25,m_1=10,m_2=0`$，曲线相同：

```math
C_G(\rho)=\frac29\rho^2+\frac29\rho^3.
```

前者连通，后者不连通。所以本命题不声称恢复图同构类型或连通性。它是已知精确函数时的代数可识别性结论；这里也没有声称有限次随机实验就能无误重建系数。

## 9. 附加结果：另一种观测也能反演

对任意非空有限简单 $`d`$-正则图（这里允许三角形及 $`d=0`$），定义
$`X_G(t)=\sum_v\mathbf1\{U_v\le t\}I_v(U)`$，
$`V_G(t)=\mathrm{Var}(X_G(t))`$。

令 $`s=1-t`$，$`p=(1-s^{d+1})/(d+1)`$。同样按中心标签积分可得

```math
V_G(t)=np-n^2p^2+
4\sum_{c=0}^d\frac{m_c}{d-c+1}
\left[p-\frac{1-s^{2d+2-c}}{2d+2-c}\right].
\tag{16}
```

证明细节：单点入选概率为 $`p`$；相邻顶点无法同时入选；非邻顶点的联合概率为
$`2\int_0^t\int_0^y(1-x)^{d-c}(1-y)^d dxdy`$，等于 (16) 方括号的 $`2/(d-c+1)`$ 倍。再展开二阶矩即得 (16)。

对 $`1\le c\le d`$，次数 $`2d+2-c`$ 不与常数、$`d+1`$、$`2d+2`$ 重叠，故

```math
m_c=\frac{(d-c+1)(2d+2-c)}4
[s^{2d+2-c}]V_G(1-s).
\tag{17}
```

$`m_0`$ 用总非邻对数减去其余 $`m_c`$ 得到。又因
$`p=t-dt^2/2+O(t^3)`$，每个非邻对的联合概率为 $`t^2+O(t^3)`$，有

```math
V_G(t)=nt-\frac{n(2+3d)}2t^2+O(t^3),
```

从而

```math
n=V_G'(0),\qquad d=-\frac{V_G''(0)/n+2}{3}.
```

所以这一阈值方差曲线也能独立恢复 $`(n,d,m_0,\ldots,m_d)`$。这里同样不将旧有峰值/峰值标签模型本身宣称为新概念。

## 10. 验证方式与实际结果

验证程序使用 Python 标准库和 `fractions.Fraction`；没有浮点容差，也没有蒙特卡洛采样。

1. **阈值公式：30 个测试图。** 对每个顶点子集，用动态规划计数该子集的全部排列前缀及局部极小点数量的一、二阶矩，再与标签落在阈值内的二项分布混合。所得整条多项式与 (16) 完全一致。测试包括环、完全图、完全二部图、非连通图、立方体、三棱柱、Petersen 图、Heawood 图，最大 14 个顶点。
2. **噪声公式：10 个测试图。** 枚举保留标签的子集。每一对局部极小事件转成有限偏序，精确计数其全部线性扩张，以获得联合概率；再展开保留概率多项式。所得整条多项式与 (2) 完全一致，并成功恢复全部 $`m_j`$。测试为 $`K_{1,1},C_4,C_5,C_6,C_7,C_8,K_{3,3},K_{4,4},2C_4`$ 和立方体。
3. **第二份独立实现：9 个测试图。** 另一审查分支独立实现了偏序计数，也全部通过；其中额外包含 Petersen 图。
4. **独立代数审查。** 两个分支分别核对了阈值反演和噪声系数，其中噪声审查使用 Hoeffding 分解，与本文条件概率证明不同。初稿阈值系数曾有一个符号笔误，审查纠正为 (17) 的正号；最终公式已通过精确检验。

以上是有限实例核验和模型间审查，不是人工同行评审，也不是 Lean/Coq 形式化认证。对任意图规模的数学依据是第 3—6 节的证明。

复现实验：在仓库根目录运行 `python verify_threshold.py`、`python verify_noise.py` 和 `python independent_noise_check.py`。前两项将结果写到 仓库根目录中的结果文件。

## 11. 文献与工具出处

- Baldi, P., Rinott, Y., Stein, C. (1989). *On the Distribution of the Number of Local Minima of a Random Function on a Graph*. [NeurIPS 官方全文](https://papers.neurips.cc/paper_files/paper/1989/file/d947bf06a885db0d477d707121934ff8-Paper.pdf)。已阅读这篇 6 页论文：相同随机标号模型，正态逼近及特定图的无阈值方差；不将这些内容作为本次原创。
- Baldi, P., Rinott, Y., Stein, C. *A Normal Approximation for the Number of Local Maxima of a Random Function on a Graph*. [Stanford 1988 年技术报告记录](https://statistics.stanford.edu/technical-reports/normal-approximation-number-local-maxima-random-function-graph)，1989 年书章 [DOI](https://doi.org/10.1016/B978-0-12-058470-3.50012-4)。确认了记录及摘要，但本次未完成长文版本逐页审读，属于查重缺口。
- Diaz-Lopez, A., Everham, L., Harris, P. E., Insko, E., Marcantonio, V., Omar, M. (2017). *Peaks on Graphs*. [arXiv:1708.08493](https://arxiv.org/abs/1708.08493)。指定峰顶点集合的标号计数。此处对其依据限于摘要和所查内容，不声称排除其全部后续文献。
- Bozeman, C., Cheng, C., Harris, P. E., Lasinis, S., Walker, S. (2024). *The Pinnacle Sets of a Graph*. [arXiv:2406.19562](https://arxiv.org/abs/2406.19562)。局部峰值的标签集合及其可实现性；独立审查分支检索了全文，未发现本处的反演陈述。
- [math1um/objects-invariants-properties](https://github.com/math1um/objects-invariants-properties)：按用户要求先检查的现成图论猜想工具库。本次没有将库内猜想冒充新命题；精确实验所需的小型偏序计数器另行实现。
- Kassis, T., Agarwal, V., He, Y., Patel, D., Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. [arXiv:2609.00065](https://arxiv.org/abs/2609.00065)。使用其中 scientific-brainstorming、paper-lookup 的工作流程；该引用只说明工具来源，不为数学结论或原创性背书。

最终状态：证明与计算已完成；绝对历史原创性未能证实。不能把这份研究笔记宣传为已经获得首创认证的数学发现。
