# GMXPolymer by Jianchuan Liu  XHU

### GMXPolymer使用教程：
https://www.simuhome.cn/forum.php?mod=viewthread&tid=110&extra=page%3D1

# 交联处理工具集

## 1. delete_itp_atom.py - ITP文件原子删除工具

### 功能说明

该工具用于从GROMACS ITP（Include Topology）文件中删除指定的原子，并自动重新编号所有相关的拓扑参数。

### 主要功能

- 解析ITP文件的各个节（[ atoms ]、[ bonds ]、[ pairs ]、[ angles ]、[ dihedrals ]）
- 删除指定的原子
- 自动删除与被删原子相关的键、角度、二面角等参数
- 重新编号所有保留的原子和相关参数
- 将被删原子的电荷平均分配到剩余原子上
- 输出修改后的ITP文件

### 使用方法

#### 程序函数

```python
from delete_itp_atom import parse_itp_file, delete_atoms_and_renumber, write_itp_file

# 1. 解析ITP文件
itp_data = parse_itp_file("input.itp")

# 2. 指定要删除的原子ID列表
atoms_to_delete = [11, 298]

# 3. 计算剩余原子数
atoms = itp_data.get('atoms')
atom_num = len(atoms) - len(atoms_to_delete)

# 4. 删除原子并重新编号
modified_data, id_mapping = delete_atoms_and_renumber(itp_data, atoms_to_delete, atom_num)

# 5. 写入新的ITP文件
write_itp_file(modified_data, "output.itp")
```

#### 直接运行脚本

修改脚本末尾的参数后直接运行：

```python
# 设置参数
itp_file = "LMW.itp"           # 输入ITP文件
atoms_to_delete = [11, 298]    # 要删除的原子ID

# 运行
python delete_itp_atom.py
```

输出文件将保存为 `modified.itp`

### 参数说明

| 参数 | 说明 |
|------|------|
| `itp_file` | 输入的ITP文件路径 |
| `atoms_to_delete` | 要删除的原子ID列表（从1开始） |

### 支持的ITP节

- [ moleculetype ] - 分子类型信息
- [ atoms ] - 原子定义
- [ bonds ] - 键参数
- [ pairs ] - 非键对参数
- [ angles ] - 键角参数
- [ dihedrals ] - 二面角参数

### 注意事项

1. 原子ID从1开始编号
2. 删除原子后，其电荷会平均分配到所有剩余原子上
3. 所有包含被删原子的键、角度、二面角都会被删除
4. 原子ID会自动重新编号，保持连续性

## 2. Breakpoint_extraction.py - 模拟断点文件提取工具

### 功能说明

该工具用于从GMXPolymer模拟的断点中提取可继续模拟的文件。当模拟中断或需要从某个检查点继续时，使用此工具可以快速整理所需的文件。

### 主要功能

- 交互式输入断点编号
- 从指定的断点目录中提取文件
- 自动创建输出目录
- 复制并重命名gro文件（重命名为 init.gro）
- 复制并重命名top文件（重命名为 top.top）
- 复制所有itp文件

### 使用方法

#### 运行脚本

```bash
python Breakpoint_extraction.py
```

#### 交互式输入

运行后会提示输入断点编号：

```
输入你想从哪一个断点提取可模拟的文件(整数): 5
```

#### 输出

脚本会创建 `BreakMD_X` 目录（X为断点编号），包含：

```
BreakMD-5/
├── init.gro    # 从 BondSteep-5/md5.gro 复制
├── top.top     # 从 BondSteep-4/B4.top 复制
└── *.itp       # 从 BondSteep-4/ 复制的所有itp文件
```

### 输出信息

运行成功后会显示：

```
!!!The all file is in the 'BreakMD_5' folder!!!
!!!The box file will named 'init.gro'; the top file will named 'top.top'!!!
```

### 文件命名规则

| 原始文件 | 目标文件 | 说明 |
|----------|----------|------|
| `BondSteep-X/mdX.gro` | `BreakMD-X/init.gro` | Gro文件重命名为 init.gro |
| `BondSteep-(X-1)/B(X-1).top` | `BreakMD-X/top.top` | Top文件重命名为 top.top |
| `BondSteep-(X-1)/*.itp` | `BreakMD-X/*.itp` | ITP文件保持原名 |

### 注意事项

1. 断点编号必须是正整数
2. Gro文件来自当前断点目录（`BondSteep-X`）
3. Top和ITP文件来自前一个断点目录（`BondSteep-(X-1)`）
4. 输出目录如果已存在会被覆盖
5. 脚本使用 `cp -rf` 命令进行文件复制

# 3. remove_molecules.py 去除交联过程中未反应的单体工具

## 功能概述

`remove_molecules.py` 是一个用于GROMACS分子动力学模拟的文件处理工具，可以同时处理拓扑文件（top）和坐标文件（gro），实现分子的批量删除，确保两个文件的一致性。

### 核心功能

1. **Top文件处理**
   - 删除 `[ molecules ]` 节中指定的分子记录
   - 删除对应的 `#include` 引用行
   - 保留注释行（以分号 `;` 开头）
   - 自动识别并删除数量为0的分子

2. **Gro文件处理**
   - 删除指定分子的所有原子坐标
   - 自动重新编号原子序号
   - 自动重新编号残基序号

3. **文件管理**
   - 将处理后的文件输出到指定文件夹
   - 自动拷贝保留的itp文件
   - 不拷贝被删除分子的itp文件

## 使用方法

### 命令行语法

```bash
python remove_molecules.py [选项]
```

### 基本示例

```bash
# 显示文件信息（不执行删除）
python remove_molecules.py -f test.top -g init.gro --info

# 删除数量为0的分子
python remove_molecules.py -f test.top -g init.gro --zero-count

# 删除指定的分子
python remove_molecules.py -f test.top -g init.gro -m B1 B2 B3

# 删除指定分子并同时删除数量为0的分子
python remove_molecules.py -f test.top -g init.gro -m B1 B2 --zero-count

# 指定输出文件夹
python remove_molecules.py -f test.top -g init.gro -m B1 -o MyOutput

# 不重新编号原子和残基
python remove_molecules.py -f test.top -g init.gro -m B1 --no-renumber
```

## 命令行参数

| 参数 | 短格式 | 默认值 | 说明 |
|------|--------|--------|------|
| `--file` | `-f` | `test.top` | 输入的top文件路径 |
| `--gro` | `-g` | `init.gro` | 输入的gro文件路径 |
| `--molecules` | `-m` | - | 要删除的分子名称列表（空格分隔） |
| `--zero-count` | - | `False` | 同时删除数量为0的分子 |
| `--output` | `-o` | `Remove` | 输出文件夹名称 |
| `--info` | - | `False` | 仅显示文件信息，不执行删除 |
| `--no-renumber` | - | `False` | 不重新编号原子和残基 |

## 输入文件要求

### Top文件格式

```
[ defaults ]
...

[ atomtypes ]
...

#include "molecule1.itp"
#include "molecule2.itp"
...

[ system ]
System Name

[ molecules ]
; molname    count
molecule1      10
molecule2       5
molecule3       0
```

### Gro文件格式

标准GROMACS gro文件格式：
- 第1行：标题
- 第2行：原子总数
- 第3行起：原子坐标（每行一个原子）
- 最后一行：盒子尺寸

## 输出文件

### 输出文件夹结构

```
Remove/
├── test.top          # 处理后的top文件
├── init.gro          # 处理后的gro文件
├── molecule1.itp     # 保留的itp文件
├── molecule2.itp     # 保留的itp文件
└── ...
```

### 处理报告

脚本运行后会输出详细的处理报告：

```
要删除的分子: ['B1', 'B2', 'B3']

处理gro文件...
  B1: 残基编号 [101]
  B2: 残基编号 [102]
  B3: 残基编号 [103]
  删除了 87 个原子
  已重新编号原子和残基

已生成新的top文件: Remove/test.top
已生成新的gro文件: Remove/init.gro

拷贝保留的itp文件到 Remove:
  已拷贝: PIP.itp -> Remove/PIP.itp
  已拷贝: TMC.itp -> Remove/TMC.itp
  ...

============================================================
处理完成!
============================================================
删除的分子: ['B1', 'B2', 'B3']
删除的原子数: 87
保留的原子数: 1623
保留的itp文件数: 8
输出文件夹: Remove
```

## 使用流程

### 典型工作流程

```
1. 查看文件信息
   └── python remove_molecules.py -f top.top -g init.gro --info

2. 确定要删除的分子
   └── 根据 --info 输出，确定删除策略

3. 执行删除操作
   └── python remove_molecules.py -f top.top -g init.gro -m B1 B2 --zero-count

4. 检查输出文件
   └── 查看 Remove/ 文件夹中的文件

5. 使用处理后的文件进行后续模拟
   └── cd Remove && gmx grompp ...
```

## 应用场景

### 场景1: 清理未反应的单体

在交联模拟中，可能存在大量未参与反应的单体分子：

```bash
# 删除所有数量为0的分子（未交联的位点）
python remove_molecules.py -f top.top -g init.gro --zero-count
```

### 场景2: 删除特定的副产物

```bash
# 删除指定的副产物分子
python remove_molecules.py -f top.top -g init.gro -m Byproduct1 Byproduct2
```

### 场景3: 创建子系统

```bash
# 只保留主要组分，删除其他所有分子
python remove_molecules.py -f top.top -g init.gro -m Water Ions --zero-count
```

## 技术细节

### 原子重新编号规则

删除原子后，脚本会自动重新编号：

1. **原子序号**: 从1开始连续编号
2. **残基序号**: 保持原有顺序，从1开始连续编号

### 残基匹配

脚本通过残基名称（residue_name）匹配要删除的分子：
- Top文件中的分子名称必须与Gro文件中的残基名称一致
- 支持一个分子对应多个残基的情况

## 注意事项

1. **文件备份**: 脚本不会覆盖原始文件，但建议在运行前备份重要文件
2. **命名一致**: 确保top文件中的分子名称与gro文件中的残基名称一致
3. **itp文件**: 被删除分子的itp文件不会被拷贝到输出文件夹
4. **注释保留**: 所有以分号开头的注释行都会被保留
5. **重新编号**: 默认会重新编号原子和残基，使用 `--no-renumber` 可禁用

## 错误处理

### 常见错误

| 错误信息 | 原因 | 解决方法 |
|----------|------|----------|
| `错误: 文件不存在 xxx` | 输入文件路径错误 | 检查文件路径是否正确 |
| `没有指定要删除的分子` | 未指定删除选项 | 使用 `-m` 或 `--zero-count` |
| `警告: 文件不存在 xxx.itp` | itp文件路径错误 | 检查itp文件是否存在 |


## 系统要求

- Python 3.6+
- GROMACS（用于后续模拟）
- 标准Python库：re, os, collections

## 作者信息

-  Jianchuan Liu, Xihua University, liujianchuan@xhu.edu.cn

## 许可

本工具集供学术研究使用。

