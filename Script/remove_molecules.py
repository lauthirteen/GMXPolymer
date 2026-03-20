#!/usr/bin/env python3
"""
GROMACS top文件和gro文件分子删除工具
功能：
1. 删除[ molecules ]中指定分子或数量为0的分子
2. 删除对应的#include行
3. 同时删除gro文件中对应的原子行
4. 保留注释行
5. 生成新的top文件和gro文件到Remove文件夹
6. 拷贝保留的itp文件到Remove文件夹

使用示例:
  # 删除数量为0的分子
  python remove_molecules.py -f test.top -g init.gro --zero-count
  
  # 删除数量为0的分子和指定的分子
  python remove_molecules.py -f test.top -g init.gro --zero-count -m PIP TMC
  # 删除指定的分子
  python remove_molecules.py -f test.top -g init.gro -m B1 B2 B3
  
  # 仅显示文件信息
  python remove_molecules.py -f test.top -g init.gro --info
"""

import os
import re
import shutil
import argparse
from typing import List, Dict, Set, Tuple


class GroFile:
    """GROMACS .gro文件解析器"""
    
    def __init__(self, filename: str):
        self.filename = filename
        self.title = ""
        self.num_atoms = 0
        self.atoms = []
        self.box_line = ""
        self._parse()
    
    def _parse(self):
        """解析gro文件"""
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        
        self.title = lines[0].strip()
        self.num_atoms = int(lines[1].strip())
        
        for i in range(2, 2 + self.num_atoms):
            line = lines[i]
            atom = {
                'residue_num': int(line[0:5]),
                'residue_name': line[5:10].strip(),
                'atom_name': line[10:15].strip(),
                'atom_num': int(line[15:20]),
                'x': float(line[20:28]),
                'y': float(line[28:36]),
                'z': float(line[36:44]),
                'has_velocity': len(line.strip()) > 44,
                'vx': 0.0,
                'vy': 0.0,
                'vz': 0.0
            }
            
            if atom['has_velocity']:
                try:
                    atom['vx'] = float(line[44:52])
                    atom['vy'] = float(line[52:60])
                    atom['vz'] = float(line[60:68])
                except:
                    atom['has_velocity'] = False
            
            self.atoms.append(atom)
        
        # box尺寸行
        if len(lines) > 2 + self.num_atoms:
            self.box_line = lines[2 + self.num_atoms]
    
    def get_residue_names(self) -> Set[str]:
        """获取所有残基名称"""
        return set(atom['residue_name'] for atom in self.atoms)
    
    def get_residues_by_name(self, mol_name: str) -> List[int]:
        """获取指定分子名称的所有残基编号"""
        residue_nums = set()
        for atom in self.atoms:
            if atom['residue_name'] == mol_name:
                residue_nums.add(atom['residue_num'])
        return sorted(list(residue_nums))
    
    def remove_residues(self, residue_nums: Set[int]) -> List[Dict]:
        """
        删除指定残基编号的原子
        
        Args:
            residue_nums: 要删除的残基编号集合
        
        Returns:
            保留的原子列表
        """
        kept_atoms = []
        for atom in self.atoms:
            if atom['residue_num'] not in residue_nums:
                kept_atoms.append(atom)
        return kept_atoms
    
    def write(self, filename: str, atoms: List[Dict] = None):
        """
        写入gro文件
        
        Args:
            filename: 输出文件名
            atoms: 原子列表（如果为None，使用原始原子列表）
        """
        if atoms is None:
            atoms = self.atoms
        
        with open(filename, 'w') as f:
            # 标题
            f.write(f"{self.title}\n")
            
            # 原子数
            f.write(f"{len(atoms):5d}\n")
            
            # 原子信息
            for atom in atoms:
                if atom['has_velocity']:
                    line = (f"{atom['residue_num']:5d}"
                           f"{atom['residue_name']:5s}"
                           f"{atom['atom_name']:5s}"
                           f"{atom['atom_num']:5d}"
                           f"{atom['x']:8.3f}"
                           f"{atom['y']:8.3f}"
                           f"{atom['z']:8.3f}"
                           f"{atom['vx']:8.4f}"
                           f"{atom['vy']:8.4f}"
                           f"{atom['vz']:8.4f}\n")
                else:
                    line = (f"{atom['residue_num']:5d}"
                           f"{atom['residue_name']:5s}"
                           f"{atom['atom_name']:5s}"
                           f"{atom['atom_num']:5d}"
                           f"{atom['x']:8.3f}"
                           f"{atom['y']:8.3f}"
                           f"{atom['z']:8.3f}\n")
                f.write(line)
            
            # box尺寸
            f.write(self.box_line)


def renumber_atoms(atoms: List[Dict]) -> List[Dict]:
    """重新编号原子"""
    for i, atom in enumerate(atoms):
        atom['atom_num'] = i + 1
    return atoms


def renumber_residues(atoms: List[Dict]) -> List[Dict]:
    """重新编号残基"""
    if not atoms:
        return atoms
    
    # 建立旧残基编号到新编号的映射
    residue_map = {}
    new_res_num = 1
    
    for atom in atoms:
        old_res_num = atom['residue_num']
        if old_res_num not in residue_map:
            residue_map[old_res_num] = new_res_num
            new_res_num += 1
        atom['residue_num'] = residue_map[old_res_num]
    
    return atoms


class TopFileProcessor:
    """top文件处理器"""
    
    def __init__(self, top_file: str):
        """
        初始化处理器
        
        Args:
            top_file: top文件路径
        """
        self.top_file = top_file
        self.lines = []
        self.includes = []
        self.molecules = []
        self.includes_section_start = -1
        self.includes_section_end = -1
        self.molecules_section_start = -1
        self.molecules_section_end = -1
        
        self._parse()
    
    def _parse(self):
        """解析top文件"""
        with open(self.top_file, 'r') as f:
            self.lines = f.readlines()
        
        in_molecules_section = False
        
        for i, line in enumerate(self.lines):
            stripped = line.strip()
            
            if stripped.startswith('#include'):
                self.includes.append({
                    'line_num': i,
                    'content': line,
                    'file': self._extract_include_file(stripped)
                })
                if self.includes_section_start == -1:
                    self.includes_section_start = i
                self.includes_section_end = i
            
            if stripped == '[ molecules ]':
                in_molecules_section = True
                self.molecules_section_start = i
                continue
            
            if in_molecules_section:
                if stripped == '' or stripped.startswith(';'):
                    continue
                
                if stripped.startswith('['):
                    in_molecules_section = False
                    self.molecules_section_end = i
                    continue
                
                parts = stripped.split()
                if len(parts) >= 2:
                    mol_name = parts[0]
                    try:
                        mol_count = int(parts[1])
                    except ValueError:
                        continue
                    
                    self.molecules.append({
                        'line_num': i,
                        'name': mol_name,
                        'count': mol_count,
                        'content': line
                    })
        
        if in_molecules_section:
            self.molecules_section_end = len(self.lines)
    
    def _extract_include_file(self, line: str) -> str:
        """从#include行提取文件名"""
        match = re.search(r'#include\s+["<](.+?)[">]', line)
        if match:
            return match.group(1)
        return ''
    
    def get_zero_count_molecules(self) -> List[str]:
        """获取数量为0的分子名称"""
        return [mol['name'] for mol in self.molecules if mol['count'] == 0]
    
    def get_all_molecules(self) -> List[str]:
        """获取所有分子名称"""
        return [mol['name'] for mol in self.molecules]
    
    def get_included_files(self) -> List[str]:
        """获取所有#include的文件"""
        return [inc['file'] for inc in self.includes if inc['file']]
    
    def remove_molecules(self, mol_names: List[str]) -> Tuple[str, Set[str]]:
        """
        删除指定的分子
        
        Args:
            mol_names: 要删除的分子名称列表
        
        Returns:
            (新的top文件内容, 保留的itp文件集合)
        """
        mol_names_set = set(mol_names)
        
        lines_to_remove = set()
        
        for mol in self.molecules:
            if mol['name'] in mol_names_set:
                lines_to_remove.add(mol['line_num'])
        
        itp_to_remove = set()
        for inc in self.includes:
            file_name = inc['file']
            mol_name = self._get_mol_name_from_file(file_name)
            
            if mol_name in mol_names_set:
                lines_to_remove.add(inc['line_num'])
                itp_to_remove.add(file_name)
        
        kept_itp_files = set()
        for inc in self.includes:
            if inc['line_num'] not in lines_to_remove and inc['file']:
                kept_itp_files.add(inc['file'])
        
        new_lines = []
        for i, line in enumerate(self.lines):
            if i not in lines_to_remove:
                new_lines.append(line)
        
        return ''.join(new_lines), kept_itp_files
    
    def _get_mol_name_from_file(self, file_path: str) -> str:
        """从文件路径提取分子名"""
        file_name = os.path.basename(file_path)
        mol_name = os.path.splitext(file_name)[0]
        mol_name = mol_name.lstrip('./')
        return mol_name
    
    def print_info(self):
        """打印文件信息"""
        print("=" * 60)
        print("Top文件信息")
        print("=" * 60)
        
        print("\n#include 文件:")
        for inc in self.includes:
            print(f"  {inc['file']}")
        
        print("\n[ molecules ] 分子信息:")
        for mol in self.molecules:
            status = " (数量为0)" if mol['count'] == 0 else ""
            print(f"  {mol['name']:15s} {mol['count']:5d}{status}")
        
        print()


def create_remove_folder(folder: str = 'Remove') -> str:
    """创建Remove文件夹"""
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"已创建文件夹: {folder}")
    return folder


def copy_itp_files(itp_files: Set[str], dest_folder: str) -> Tuple[int, int]:
    """拷贝itp文件到目标文件夹"""
    copied = 0
    failed = 0
    
    for itp_file in itp_files:
        itp_file = itp_file.lstrip('./')
        
        if os.path.exists(itp_file):
            dest_file = os.path.join(dest_folder, os.path.basename(itp_file))
            shutil.copy2(itp_file, dest_file)
            print(f"  已拷贝: {itp_file} -> {dest_file}")
            copied += 1
        else:
            print(f"  警告: 文件不存在 {itp_file}")
            failed += 1
    
    return copied, failed


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='GROMACS top文件和gro文件分子删除工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 删除数量为0的分子
  python remove_molecules.py -f test.top -g init.gro --zero-count
  
  # 删除指定的分子
  python remove_molecules.py -f test.top -g init.gro -m B1 B2 B3
  
  # 仅显示文件信息
  python remove_molecules.py -f test.top -g init.gro --info
        """
    )
    
    parser.add_argument('-f', '--file', default='test.top',
                        help='输入的top文件 (默认: test.top)')
    parser.add_argument('-g', '--gro', default='init.gro',
                        help='输入的gro文件 (默认: init.gro)')
    parser.add_argument('-m', '--molecules', nargs='+', default=[],
                        help='要删除的分子名称列表')
    parser.add_argument('--zero-count', action='store_true',
                        help='同时删除数量为0的分子')
    parser.add_argument('-o', '--output', default='Remove',
                        help='输出文件夹 (默认: Remove)')
    parser.add_argument('--info', action='store_true',
                        help='仅显示文件信息，不执行删除')
    parser.add_argument('--no-renumber', action='store_true',
                        help='不重新编号原子和残基')
    
    args = parser.parse_args()
    
    # 检查文件是否存在
    if not os.path.exists(args.file):
        print(f"错误: 文件不存在 {args.file}")
        return
    
    if not os.path.exists(args.gro):
        print(f"错误: 文件不存在 {args.gro}")
        return
    
    # 创建处理器
    processor = TopFileProcessor(args.file)
    
    # 如果只是显示信息
    if args.info:
        processor.print_info()
        
        # 显示gro文件信息
        gro = GroFile(args.gro)
        print("=" * 60)
        print("Gro文件信息")
        print("=" * 60)
        print(f"文件: {args.gro}")
        print(f"原子数: {gro.num_atoms}")
        print(f"残基名称: {gro.get_residue_names()}")
        print()
        return
    
    # 确定要删除的分子
    mol_to_remove = set(args.molecules)
    
    if args.zero_count:
        zero_mols = processor.get_zero_count_molecules()
        mol_to_remove.update(zero_mols)
        print(f"数量为0的分子: {zero_mols}")
    
    if not mol_to_remove:
        print("没有指定要删除的分子。使用 -m 指定分子名称或 --zero-count 删除数量为0的分子。")
        processor.print_info()
        return
    
    print(f"\n要删除的分子: {list(mol_to_remove)}")
    
    # 处理top文件
    new_top_content, kept_itp_files = processor.remove_molecules(list(mol_to_remove))
    
    # 处理gro文件
    print("\n处理gro文件...")
    gro = GroFile(args.gro)
    
    # 找到要删除的残基编号
    residues_to_remove = set()
    for mol_name in mol_to_remove:
        residue_nums = gro.get_residues_by_name(mol_name)
        residues_to_remove.update(residue_nums)
        print(f"  {mol_name}: 残基编号 {residue_nums}")
    
    # 删除原子
    kept_atoms = gro.remove_residues(residues_to_remove)
    print(f"  删除了 {gro.num_atoms - len(kept_atoms)} 个原子")
    
    # 重新编号
    if not args.no_renumber:
        kept_atoms = renumber_atoms(kept_atoms)
        kept_atoms = renumber_residues(kept_atoms)
        print("  已重新编号原子和残基")
    
    # 创建输出文件夹
    output_folder = create_remove_folder(args.output)
    
    # 保存新的top文件
    output_top = os.path.join(output_folder, os.path.basename(args.file))
    with open(output_top, 'w') as f:
        f.write(new_top_content)
    print(f"\n已生成新的top文件: {output_top}")
    
    # 保存新的gro文件
    output_gro = os.path.join(output_folder, os.path.basename(args.gro))
    gro.write(output_gro, kept_atoms)
    print(f"已生成新的gro文件: {output_gro}")
    
    # 拷贝保留的itp文件
    print(f"\n拷贝保留的itp文件到 {output_folder}:")
    copied, failed = copy_itp_files(kept_itp_files, output_folder)
    
    # 打印总结
    print("\n" + "=" * 60)
    print("处理完成!")
    print("=" * 60)
    print(f"删除的分子: {list(mol_to_remove)}")
    print(f"删除的原子数: {gro.num_atoms - len(kept_atoms)}")
    print(f"保留的原子数: {len(kept_atoms)}")
    print(f"保留的itp文件数: {copied}")
    if failed > 0:
        print(f"未找到的itp文件数: {failed}")
    print(f"输出文件夹: {output_folder}")


if __name__ == '__main__':
    main()
