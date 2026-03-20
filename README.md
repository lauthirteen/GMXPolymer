# GMXPolymer by Jianchuan Liu  XHU  2024-03-26

使用教程：
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



## 系统要求

- Python 3.6+
- GROMACS（用于后续模拟）
- 标准Python库：re, os, collections

## 作者信息

-  Jianchuan Liu, Xihua University, liujianchuan@xhu.edu.cn

## 许可

本工具集供学术研究使用。

