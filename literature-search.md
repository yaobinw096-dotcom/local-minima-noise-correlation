# 原创性定向检索记录

日期：2026-10-10。检索对象：正则图随机局部极小点数量的重采样协方差多项式，以及该曲线对非邻顶点对共同邻居数量分布的反演。

结论限于：本次返回并检查的文献中未发现同一反演陈述。**未证明“此前不存在”，也未完成系统文献综述。** 不能从检索无命中推出首次提出。

## 1. GitHub：先检查现有工作与工具

初始查询：

```text
site:github.com combinatorics conjecture graph invariants experimental mathematics
site:github.com cyclic words combinatorics enumeration conjecture
```

找到并阅读了 [math1um/objects-invariants-properties](https://github.com/math1um/objects-invariants-properties) 的 README。该项目收集图、图不变量、性质、已知定理及自动猜想工具。本次未使用其中已有猜想充当原创结果。

后续定向查询：

```text
site:github.com "local minima" "noise" "codegree"
```

返回结果中未发现该精确命题。GitHub 的索引覆盖有限，不等于逐库代码检索或全文数据库检索。

## 2. 网页定向检索

代表性实际查询如下；搜索引擎对部分过窄查询返回了无关近似结果，未将其计作“同题零篇文献”的证据。

```text
random local minima graph variance regular triangle free independent set
"random peaks" graph variance regular
"local minima" "variance" "common neighbors" graph
"Caro" "Wei" random independent set variance
"local minima" "regular graph" "variance"
"local maxima" "variance" "graph" random labeling
"Boppana" "variance"
"random labelling" "peaks" graph variance
"local maxima" graph "threshold" "variance"
"local minima" "codegree"
"local maxima" graph "reconstruct" variance
"local maxima" "noise stability" graph
"local minima" "resampling" graph covariance
"random function on a graph" covariance
"local extrema" "noise sensitivity" graph
"local maxima" "Hoeffding" graph
"local minima" "common neighbours" variance
"local minima" "covariance" "resampling" graph vertices
"peaks on graphs" noise covariance
"local minima" "binomial" "Rinott"
```

采用英文表达、local minima/local maxima 对偶表达及 common neighbors/common neighbours/codegree 等变体，避免只检索本次命名。

## 3. arXiv API 的可复现查询

端点：`https://export.arxiv.org/api/query`。

参数：

```text
search_query=(all:"local minima" OR all:"local maxima") AND (all:graph) AND (all:noise OR all:variance OR all:threshold)
start=0（默认）
max_results=20
```

实际响应：`total_results=18`，`returned=18`，并核对了服务端回显查询。使用 paper-lookup 附带 `arxiv_atom.py` 解析原始 Atom XML，未把 HTTP 成功状态当作检索成功的充分条件。

18 条主要涉及优化、SLAM、Ising 系统等同名术语；标题与摘要筛查未见本命题。该查询甚至没有覆盖到最接近的 1989 年经典文献，因此其未命中只提供很弱的新颖性证据。压缩包附 `arxiv_search.json`，包含实际返回的题录和摘要。

## 4. 找到的最近似文献

| 文献 | 本次检查范围 | 与候选的关系 |
|---|---|---|
| [Baldi–Rinott–Stein，1989，On the Distribution of the Number of Local Minima of a Random Function on a Graph](https://papers.neurips.cc/paper_files/paper/1989/file/d947bf06a885db0d477d707121934ff8-Paper.pdf) | 官方 6 页全文 | 相同随机标号模型、正态逼近、部分图族的方差。未在该文中看到整条重采样协方差曲线的共同邻居分布反演。 |
| [A Normal Approximation for the Number of Local Maxima of a Random Function on a Graph](https://statistics.stanford.edu/technical-reports/normal-approximation-number-local-maxima-random-function-graph) | Stanford 1988 技术报告题录与 1989 书章摘要；定位到 PDF，但未完成长文逐页审读 | 重要查重缺口；不能宣称已彻底排除该长文及其后续工作的等价结果。 |
| [Peaks on Graphs，2017](https://arxiv.org/abs/1708.08493) | 摘要、相关出版记录 | 指定峰顶点集合的标号计数；说明峰统计是已有研究对象。 |
| [The Pinnacle Sets of a Graph，2024](https://arxiv.org/abs/2406.19562) | 摘要；独立审查分支检索 HTML 全文 | 峰值标签集合、可实现性、偏序；未找到本处的方差/噪声反演陈述。 |
| [Variance of random greedy independent sets in triangle-free graphs，2026](https://arxiv.org/abs/2609.14826) | 搜索返回的论文摘要 | 使用迭代 random greedy independent set，随机变量与本命题的一轮局部极小点计数不同；标题相近不能视为同一结果。 |

## 5. 明确没有完成的工作

没有完整检索 MathSciNet、zbMATH、全部学位论文、全部非英文材料、所有引用链或未公开手稿；没有获得人类领域专家审查。Stanford 长报告的官方 PDF 下载尝试在 25 秒后超时、收到 0 字节，因此不能把它记为已读全文。还需警惕结果以噪声敏感性、图上峰值、随机能量模型、Matérn 型筛选或其他等价表述出现。

本次候选使用标准条件概率、二项式恒等式与三角线性反演，证明方法本身不宣称原创。数学证明说明陈述成立；它无法证明历史优先权。若后续发现相同或更一般的已知命题，应将当前结果降格为独立重发现或已知结果的推论。
