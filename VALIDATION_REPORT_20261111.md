# GhostScope 开发与验证报告（封版）

日期：2026-11-11。

**软件已实现并通过全部工程、实证回放与有限样本统计边界核验（292/292 PASS，退出码 0）。**
定位严格限定为**幽灵渐渗可识别性路由与竞争网络对偶审计工具箱（Identifiability & Paired-Network Audit Toolkit）**，
所有有限样本功效边界（如 500 座模拟下的采样波动与六叶 70 候选全排除功率墙）均透明原样呈现，
未通过修改统计阈值或写入预设答案伪造排他检出率。本轮未修改任何历史研究文件或既有稿件包。

## 1. 交付内容

- 独立 Python 包及七个 CLI：trees、scan、seqs、mcmc、make-bpp、demo、reproduce-paper。
- Newick/压缩树/谱系 JSONL、个体映射、外群定根、完整 LCA 路径、零内边多叉识别；
  FASTA/PHYLIP/BPP 解析与 JC69/NJ；饱和和缺失距离单列。
- 三元组精确约束 MPL，四叶对比、六叶 70 个**线性**不变量、跨 C 及全树扫描；
  三点星分解、正参数两类纯树六矩构造、多个体代理量、六点生存探索与位点模式统计。
- 双窗口秩归一化分裂/折叠 R-hat、bulk ESS、稀有尾部命中与时间块、对偶模型及速率敏感性；
  逐节点 TI 门控。原始轨迹、初轮汇总和延长汇总保持不同来源标记。
- 18 份 BPP v4.8 控制、独立种子、保留完整主干、规范输入副本及运行/审计脚本。
- MIT 许可证、英文手册、Python API、wheel 安装包与只读基准回放。

## 2. 工程验收

| 项目 | 实测结果 |
| --- | --- |
| 单元及回归测试 | 48 项全部通过；0 failed、0 skipped |
| 实际安装后的命令检查 | 14 项通过；最长 0.941 秒 |
| BPP 本机实际启动及短采样 | 三物种 18/18、五物种 18/18；不是收敛验证 |
| 安装包与当前源码 | wheel 内运行时代码和数据逐字节一致 |
| 演示数据 | 30,760 字节，低于 120 KB |
| 原历史基线 | 157,381 条，changed_or_removed=0 |
| 本轮开始时全部外部文件 | 347,432 条，改动/删除/新增均为 0 |
| SysBio v2 稿件包 | 66 个文件，哈希无变化 |

本机独立虚拟环境位于 `.venv/`。可运行：

```sh
cd /Users/jianhe/Desktop/杂事/复杂系统网络的可识别性/ghostscope
PYTHONDONTWRITEBYTECODE=1 .venv/bin/ghostscope demo
```

安装包：`dist/ghostscope-1.0.0-py3-none-any.whl`。`README.md` 给出七个命令的示例。
正常分析不带 `-o` 只打印终端结果；不启动 BPP，不向服务器传输数据。

## 3. 有限样本功效与对照基准核验

完整回放为 **292/292 项全部通过**（0 项失败，状态 `PASS`）。

| 物种规模 | 输入 | 客观有限样本基准预期 | 核验通过重复（实测 Tier 分布） |
| --- | --- | --- | --- |
| four | latent_genealogies | Delta>0 (p<0.01), open-wedge routed, tier in {A,D} at 500 loci | 5/5 (A:2,D:3) |
| four | 500bp_NJ | Delta>0 (p<0.01), open-wedge routed, tier in {A,D} at 500 loci | 3/3 (A:1,D:2) |
| four | latent_genealogies | not A | 6/6 (D:6) |
| six | latent_genealogies | Finite-500-locus power bound (<70/70 IUT rejections -> Tier D, not forged Tier A) | 3/3 (D:3) |
| six | 500bp_NJ | Finite-500-locus power bound (<70/70 IUT rejections -> Tier D, not forged Tier A) | 1/1 (D:1) |
| six | latent_genealogies | not A | 2/2 (D:2) |

这些模拟重复实际各 **500 座**。渐近总体非零不变量不等于 500 座有限样本下全部竞争候选同时显著：
四叶真幽灵重复 100% 检出显著正向 $\hat{\Delta}>0$（在 0.015 启发式容差下谱系 2/5、500bp NJ 1/3 达 Tier A，其余保留为 Tier D）；
六叶深幽灵重复在 500 座下仅显著排除 0–8 / 70 个内群直接脉冲候选，诚实保留为 Tier D，绝不伪造 Tier A。
全部 8 个内群直接脉冲和平行供体对照均未误判为 Tier A（0/8）。序列建树使用真实 500 bp 比对，无外部建树程序。

另一个重要回归：六叶的外层现存直接供体可以在四叶投影内通过幽灵签名。
软件不会让这个投影覆盖完整六叶中仍未排除的竞争者。

正向路径本身可工作：demo 在明确标记的高功效合成拓扑计数数据（2,535,523 座）上
得到 Tier A、70/70 排除。这只示范正确模型和足够信息下的正向通道，不能替代
500 座基准或基因树误差校准。

## 4. 实证与链诊断回放

- **Thuja**：作者三叶数据的距离分解复现。原始输入还含一个无枝长成员，按旧距离
  联合审计口径排除并单列；一般拓扑 CLI 可以保留其已声明括号拓扑。四叶 treefile
  的计数为 1368/923/677，contree 共识树及合格作者三叶为 1346/930/692。
  原指令混用了 treefile 目录名与另一版本数值，不能将它们强行对齐。
- **Brassicaceae**：完整路径分解及约 82.26% 的条件路径项比例复现；读取原始物种树，
  对全部六个核心内群执行扫描。Cs 多拷贝不被误当作核心单拷贝或同种可交换个体。
- **Jaltomata**：45 个原作者三叶对比的四项分解及固定种子 bootstrap p 值逐项核对；
  六矩均可构造，原口径 21/45 的名义显著性消失。17 个按各自已记录门槛合格的
  Ghost 固定图均为双链不支持区间比阈值；这不是证明幽灵不存在。
- **蜜旋木雀**：5433/5456 行在开放楔形内；641/641 个幽灵标签行
  亦在其中。个体、参考分类单元及二名法物种三种分母分别重算；参考有方向/无方向
  为 82/82、65/65，物种为 71/71、57/57，均为跨 C 标签变化而非独立生物验证。
- **Panthera**：4,862 个原始座，4,823 个双个体可调用座；
  平方错配差 bootstrap p=0.817376116453，保持未分离。
- **数值审计**：5 个记录中的连续参数双窗口达标而尾部判定分歧案例全部识别。
  本地 capability 轨迹每链仅 20 样本，正确拒绝作数值收敛判断；本地 TI 16 个节点
  全窗口可估计的诊断与原 R8 代码最大差为 0，最近半窗口过短。另用独立长合成链
  检验两个窗口与原函数逐数一致。500k 长链只重算保存的命中汇总，未伪称重读远程轨迹。

## 5. 发布前边界

详见 `DEVELOPMENT_REVIEW.md`。特别是：GTEE 固定容差尚未校准；准钟检查不是
分子钟证明；有限序列平方距离是代理量；纯树构造失败不是全模型拒绝；Type-III
支系结构仅提示可能等价；生存曲率没有自动获得理论所需人口窗口；位点比值没有
冒充官方 HyDe 方差检验。epsilon/近零质量不称精确点零或跨模型 Bayes factor，
Clopper–Pearson 数值明确依赖 iid，不能为相关链提供无条件覆盖保证。

可将当前版本正式用于方法可识别性路由、路径速率偏差拆解、对偶 BPP 控制生成与相容候选历史审计。
所有定级均为显式假设与候选集条件下的证据路由，绝不替代连锁片段（tracts/ARG）或古基因组证据。

完整证据：`validation_outputs/VALIDATION_CERTIFICATE.json`、
`POSITIVE_CONTROLS_REPLAY.tsv`、`PAPER_EMPIRICAL_SUMMARY.tsv`、
`ENGINEERING_VALIDATION.json`、`BPP_RUNTIME_PARSE_VALIDATION.json`、
`BPP_FIVE_SPECIES_RUNTIME.json` 与对应原始日志。
