# Local minima under partial resampling

## 局部极小点计数的重采样关联：反例与精确公式

这是一个可复现的数学研究记录。**“所有无三角正则图都有同一极限关联曲线”的普适性猜想已被精确反例否定。** 本仓库保留反例、成立的有限图公式，以及共同邻居数量分布的反演证明。

基础随机标号模型与重采样操作是已有工具；这里的精确反演命题是否已有等价文献，尚未完成认证。结果由 Codex 辅助推导与审查，未经人工同行评审或证明助手形式化认证。

### 实验定义

给图的每个顶点独立抽取 $U_v\sim\mathrm{Unif}(0,1)$。令 $X_G$ 为标号小于所有邻居的顶点数量。随后，每个标号以概率 $p$ 保留，否则独立重抽，得到计数 $X_G^{(p)}$。

考察 $C_G(p)=\operatorname{Cov}(X_G,X_G^{(p)})$ 和 $\rho_G(p)=C_G(p)/C_G(1)$。本项目的主结论适用于非空有限简单无向、无三角形的 $d$-正则图，$d\ge2$。

### 反例

四环和五环的精确相关曲线为

$$
\rho_{C_4}(p)=\frac35p^2+\frac25p^3,
\qquad
\rho_{C_5}(p)=\frac12p^2+\frac12p^3.
$$

取 $G_m=5mC_4$、$H_m=4mC_5$。两者都有 $n=20m$、$d=2$，但各连通分支独立，因此相关曲线不随 $m$ 改变，且

$$
\rho_{G_m}(p)-\rho_{H_m}(p)=\frac1{10}p^2(1-p)>0
\quad(0<p<1).
$$

差异在 $n\to\infty$ 时仍保留，无法解释为 $1/n$ 修正。详见 [反例与修正命题](counterexamples.md)，其中也给出连通的三正则反例族。

### 保留下来的精确结果

对每个允许的有限图，存在依赖于图的非负权重，使

$$
\rho_G(p)=\sum_{k=2}^{d+1}w_k(G)p^k,
\qquad \sum_kw_k(G)=1,
\qquad w_2(G)>0,\ w_{d+1}(G)>0.
$$

因此小 $p$ 时有二次指数，但系数为 $w_2(G)$；$\rho_G(p)/p^2$ 严格递增，且 $0<\rho_G(p)<p^2$。完整协方差曲线还能恢复 $n,d$ 及非邻顶点对共同邻居数量的分布。详见 [完整公式和证明](exact-results.md)。

### 复现

需要 Python 3.10 或更新版本，仅使用标准库，无需安装依赖。在仓库根目录运行：

```text
python counterexamples.py
python verify_noise.py
python verify_threshold.py
python independent_noise_check.py
```

- `counterexamples.py`：独立偏序计数核验四环、五环，并检查反例和修正后的形状性质。
- `verify_noise.py`：10 个图的整条噪声多项式及共同邻居分布反演。
- `verify_threshold.py`：30 个图的阈值方差多项式。
- `independent_noise_check.py`：另一份偏序计数实现，9 个图，包括 Petersen 图。

计算使用精确有理数，不采用蒙特卡洛或浮点容差。仓库根目录中的结果文件 保存实际运行结果。有限测试用于检查实现；一般规模的依据是文档中的数学证明。

### 文献

- Baldi, Rinott, Stein (1989), [On the Distribution of the Number of Local Minima of a Random Function on a Graph](https://papers.neurips.cc/paper_files/paper/1989/file/d947bf06a885db0d477d707121934ff8-Paper.pdf)：基础模型、正态逼近及部分图族的方差。
- Ryan O'Donnell, [Analysis of Boolean Functions 课程材料](https://www.cs.cmu.edu/~odonnell/boolean-analysis/homework2.pdf)：保留坐标、否则独立重采样的标准噪声操作。
- Diaz-Lopez et al. (2017), [Peaks on Graphs](https://arxiv.org/abs/1708.08493)。
- Bozeman et al. (2024), [The Pinnacle Sets of a Graph](https://arxiv.org/abs/2406.19562)。

更详细的检索范围与缺口见 [文献检索记录](literature-search.md)。未找到相同表述不构成首创认证。
