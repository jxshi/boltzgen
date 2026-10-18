# GLP-1R 本地 BoltzGen 配置

本目录提供两个可直接用于 `peptide-anything` 的设计规格：

- `glp1r_peptide_central.yaml`：主批次，28–32 aa。
- `glp1r_peptide_wide.yaml`：扩展批次，25–35 aa；建议在主批次方向检查后再运行。

根据第一轮结构复核，新增两套下一轮配置：

- `glp1r_peptide_native_pose.yaml`：推荐优先运行。保留 6X18 中天然 GLP-1 的 H7–D15 序列和位姿，重新设计其余位置。
- `glp1r_peptide_agonist_motif.yaml`：de novo 探索方案。固定 N 端 `HAEGTFTSD`，设计后续 19–23 aa，并排除两个实验结构缺失环作为结合位点。

第一轮最终 10 个候选全部以 C 端进入天然 H7 激活口袋，因此不建议直接扩大原始 `central` 批次。下一轮应先运行 `native_pose`，并用全新的输出目录。

文件名特意不以 `_数字.yaml` 结尾；BoltzGen 将该格式保留给内部文件索引。

两者读取 `reference/6X18_R_canonical.cif`。克隆仓库后进入本目录运行即可。配置中的受体链为 R；目标结构不包含天然 GLP-1、G 蛋白或其他辅助链。

## 编号

BoltzGen的设计规格使用从1开始的mmCIF `label_seq_id`，不是PDB作者编号，也不是云端API的0基索引：

| 人源GLP-1R编号 | BoltzGen YAML编号 | 云端API编号 |
|---:|---:|---:|
| Q234 | 206 | 205 |
| W306 | 278 | 277 |
| R380 | 352 | 351 |

完整映射在 `reference/residue_mapping.tsv`。三个位置是稀疏的期望接触位点。它们不会强制每个位置同时接触，也不会强制肽 N 端进入口袋。

## 推荐执行顺序

下一轮天然位姿 smoke test：

```bash
boltzgen check glp1r_peptide_native_pose.yaml
boltzgen run glp1r_peptide_native_pose.yaml \
  --output runs/glp1r_native_pose_smoke \
  --protocol peptide-anything \
  --num_designs 100 \
  --budget 20 \
  --devices 2
```

下一轮 N 端基序引导的 de novo smoke test：

```bash
boltzgen check glp1r_peptide_agonist_motif.yaml
boltzgen run glp1r_peptide_agonist_motif.yaml \
  --output runs/glp1r_agonist_motif_smoke \
  --protocol peptide-anything \
  --num_designs 200 \
  --budget 20 \
  --devices 2
```

以下原始配置保留用于结果追溯：

先确认版本并检查规格：

```bash
git clone https://github.com/jxshi/boltzgen.git
cd boltzgen/glp1r_peptide
boltzgen --version
boltzgen check glp1r_peptide_central.yaml
```

打开`check`产生的CIF，确认显示的目标是受体R链，并且高亮位置对应Q234、W306和R380。官方文档明确建议在计算前做这一步。

小批次技术与姿态测试：

```bash
boltzgen run glp1r_peptide_central.yaml \
  --output runs/glp1r_28_32_smoke \
  --protocol peptide-anything \
  --num_designs 50 \
  --budget 10 \
  --devices 2
```

确认结果格式、GPU运行和结合方向后，再进行主批次：

```bash
boltzgen run glp1r_peptide_central.yaml \
  --output runs/glp1r_28_32_prod \
  --protocol peptide-anything \
  --num_designs 10000 \
  --budget 100 \
  --devices 2
```

当前命令按两张可用 GPU 配置为 `--devices 2`。中断后用完全相同的命令并加 `--reuse` 恢复，避免换输出目录重新开始。

可先只生成并检查解析后的流水线配置：

```bash
boltzgen configure glp1r_peptide_central.yaml \
  --output runs/glp1r_28_32_configured \
  --protocol peptide-anything \
  --num_designs 50 \
  --budget 10 \
  --devices 2

boltzgen execute runs/glp1r_28_32_configured
```

扩展长度批次：

```bash
boltzgen run glp1r_peptide_wide.yaml \
  --output runs/glp1r_25_35_prod \
  --protocol peptide-anything \
  --num_designs 5000 \
  --budget 50 \
  --devices 2
```

官方说明`--num_designs`是中间设计数量，`--budget`是经过质量和多样性筛选后的最终集合数量。这里的50/10用于流程和姿态验证，10000/100与5000/50是起步建议，不是模型要求；可按服务器吞吐量调整。

## 结果检查

主要结果在输出目录的`final_ranked_designs/`；完整指标通常在`final_ranked_designs/all_designs_metrics.csv`和`final_designs_metrics_<budget>.csv`，复合物结构来自`intermediate_designs_inverse_folded/refold_cif/`。

前期云端小批次显示，仅靠 Q234/W306/R380 引导时，10 条中有 7 条以 C 端进入天然配体 N 端口袋。因此本地 smoke test 首先应检查：

1. 肽N端还是C端靠近天然GLP-1的H7位置；
2. Q234/W306/R380的接触来自肽的哪一端；
3. 是否依赖6X18未解析的受体130–135或339–343位；
4. 与N63/N82/N115潜在糖链是否过近；
5. 预测复合物是否有明显原子碰撞及异常受体位移。

BoltzGen的排序指标用于计算筛选，不能证明实验亲和力或GLP-1R激动活性。即使方向与天然肽一致，仍需cAMP/Gs剂量反应实验区分激动剂、部分激动剂和拮抗剂。

## 当前接口依据

配置格式与命令依据BoltzGen官方[设计规格说明](https://github.com/HannesStark/boltzgen/blob/main/example/README.md)、[peptide binding-site示例](https://github.com/HannesStark/boltzgen/blob/main/example/vanilla_peptide_with_target_binding_site/beetletert.yaml)及[主README](https://github.com/HannesStark/boltzgen)。官方当前说明`peptide-anything`默认在逆折叠阶段避免生成Cys；仍应在最终序列中重新核对。

