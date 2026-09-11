# -% ===== 5.1 问题 1 的模型建立与求解 =====
\subsection{问题 1：时频冲突检测模型}

\subsubsection{时频资源的统一表示}

按题目要求，频域以 $\Delta f$ 为基本单位离散为 100 个频段，记可用频段集合为
$\mathcal{F}=\{0,1,\dots,99\}$（第 $f$ 个频段对应 $[f_0+f\Delta f,\ f_0+(f+1)\Delta f)$）；
时域以 $\Delta t$ 为基本单位，记 $\mathcal{T}=\{0,1,\dots,T-1\}$，其中
$T=\left\lceil\max_i\bigl(t_i^++(n_i-1)\tau_i\bigr)\right\rceil$ 为全部计划的最后结束时刻。
在附件 1 数据下 $T=621$。

对第 $i$ 个用频计划，定义其\textbf{时空资源占用集}

\begin{equation}
  O_i=F_i\times\mathcal{T}_i,\qquad
  F_i=\{f_i^-,f_i^-+1,\dots,f_i^+-1\},\qquad
  \mathcal{T}_i=\bigcup_{k=0}^{n_i-1}\bigl\{t_i^-+k\tau_i,\ t_i^-+k\tau_i+1,\dots,t_i^++k\tau_i-1\bigr\}.
  \label{eq:occupancy}
\end{equation}

式 \eqref{eq:occupancy} 是全文的核心表示：频段方向为一个连续区间，
时间方向为 $n_i$ 个等间隔复制的区间之并。由此定义\textbf{占用计数函数}

\begin{equation}
  \chi(f,t)=\#\{i:\ (f,t)\in O_i\},\qquad (f,t)\in\mathcal{F}\times\mathcal{T},
  \label{eq:chi}
\end{equation}

它是问题 1 的统计量与问题 3 的容量约束所共用的核心量。

\subsubsection{冲突判据与有限枚举化}

题目规定：当两个用频计划的时间区间与频段区间都存在交叠时判定为时频冲突。据此定义
\textbf{时频冲突图} $G_c=(V,E_c)$，其中 $V$ 为全部用频计划，

\begin{equation}
  \{i,j\}\in E_c\iff F_i\cap F_j\ne\varnothing\ \wedge\ \mathcal{T}_i\cap\mathcal{T}_j\ne\varnothing .
  \label{eq:edge}
\end{equation}

式 \eqref{eq:edge} 的频段部分为区间相交，可直接判定；时间部分涉及两个区间族，需要下述命题。

\begin{quote}
\textbf{命题 1（周期占用相交判据）}\quad
设 $d_i=t_i^+-t_i^-$，则 $\mathcal{T}_i\cap\mathcal{T}_j\ne\varnothing$ 当且仅当存在整数
$k_i\in[0,n_i)$、$k_j\in[0,n_j)$ 使得
\begin{equation}
  \max\bigl(t_i^-+k_i\tau_i,\ t_j^-+k_j\tau_j\bigr)
  <\min\bigl(t_i^++k_i\tau_i,\ t_j^++k_j\tau_j\bigr).
  \label{eq:criterion}
\end{equation}
\end{quote}

\textbf{证明}：$\mathcal{T}_i$ 是 $n_i$ 个左闭右开区间之并。两个并集相交当且仅当存在一对
组成区间相交（并集的元素必属于其中某个区间）。对左闭右开区间 $[a,a+d_a)$、$[b,b+d_b)$，
相交 $\iff b-a<d_a$ 且 $a-b<d_b$，等价于 $\max(a,b)<\min(a+d_a,b+d_b)$，
即式 \eqref{eq:criterion}。$\square$

命题 1 把"两个周期区间族是否相交"化为至多 $n_in_j$ 次整数比较；在附件 1 中
$n_in_j\le 12\times12=144$，故检测所需比较次数不超过
$\binom{150}{2}\times144\approx1.6\times10^6$，可在毫秒级完成，且\textbf{结果为精确判定}。

\begin{quote}
\textbf{命题 2（频段剪枝）}\quad 若 $F_i\cap F_j=\varnothing$，则 $\{i,j\}\notin E_c$。
\end{quote}

\textbf{证明}：由 $O_i\cap O_j\subseteq(F_i\cap F_j)\times(\mathcal{T}_i\cap\mathcal{T}_j)$，
频段交集为空时 $O_i\cap O_j=\varnothing$。$\square$

由命题 2，检测可分两步：先构造\textbf{频段重叠图} $G_f=(V,E_f)$，
再仅对 $E_f$ 中的对执行式 \eqref{eq:criterion} 的时间判定。附件 1 中
$|E_f|=1\,352$，时间判定次数由 $11\,175$ 降至 $1\,352$，计算量下降 88\%。

\subsubsection{检测算法}

算法 1 给出问题 1 的完整流程。

\begin{table}[!htbp]
  \centering
  \caption{问题 1 时频冲突检测算法}
  \label{tab:algo1}
  \begin{tabular}{p{0.92\linewidth}}
    \toprule[1.5pt]
    \textbf{算法 1}\quad 时频冲突检测 \\
    \midrule[1pt]
    \textbf{输入}：用频计划集合 $\{(F_i,T_i,\tau_i,n_i)\}_{i=1}^{N}$，$N=150$ \\
    \textbf{输出}：冲突对集合 $E_c$；统计量 $|E_c|$、度分布、连通分量、$\chi$ 分布 \\
    \midrule[1pt]
    1.\quad 对 $\mathcal{F}$ 上所有计划对做区间相交判定，得频段重叠图 $G_f$ \\
    2.\quad 初始化 $E_c\leftarrow\varnothing$ \\
    3.\quad \textbf{for each} $\{i,j\}\in E_f$ \textbf{do} \\
    4.\quad\quad \textbf{for} $k_i=0$ \textbf{to} $n_i-1$ \textbf{do} \\
    5.\quad\quad\quad 由式 \eqref{eq:criterion} 计算 $k_j$ 的可行区间
    $\bigl[\lfloor(s_i-t_j^-{-}d_j)/\tau_j\rfloor{+}1,\ \lceil(e_i-t_j^-)/\tau_j\rceil{-}1\bigr]
    \cap[0,n_j)$ \\
    6.\quad\quad\quad \textbf{if} 可行区间非空 \textbf{then} $E_c\leftarrow E_c\cup\{\{i,j\}\}$；
    跳出内层循环 \\
    7.\quad 由式 \eqref{eq:chi} 统计占用计数函数 $\chi$ 的分布与极值 \\
    8.\quad 用并查集求 $G_c$ 的连通分量，统计最大连通分量规模 \\
    9.\quad \textbf{return} $E_c$ 及全部统计量 \\
    \bottomrule[1.5pt]
  \end{tabular}
\end{table}

步骤 5 中的可行区间由式 \eqref{eq:criterion} 直接推出：需要
$t_j^-+k_j\tau_j<t_i^++k_i\tau_i$ 且 $t_i^-+k_i\tau_i<t_j^++k_j\tau_j+d_j$，
两式对 $k_j$ 解出上下界，再与 $[0,n_j)$ 求交。该写法把内层循环从 $n_j$ 次比较压缩为
$O(1)$ 次边界计算，是本文实现的加速要点。

\subsubsection{求解结果}

\input{5.1.1.问题1结果.tex}

\subsubsection{独立性验证}

为避免实现缺陷导致结论错误，本文用三种彼此独立的算法复算冲突总数：
精确解析枚举（算法 1）、逐时刻集合求交的暴力法（以 $\Delta t/2$ 为步长把每次占用展开为
时刻集合后求交）、以及基于 NumPy 的向量化枚举。三者在附件 1 上得到的冲突对数与冲突对
集合\textbf{完全一致}，说明检测结果不含实现误差。

