# -*- coding: utf-8 -*-
import re
import os
import argparse
from collections import defaultdict

##################################################################
# Gro文件处理函数
##################################################################

def parse_gro_file(gro_file):
    """
    解析GROMACS GRO文件
    
    返回:
    title: 标题行
    atoms: 原子信息列表
    box: box尺寸
    """
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    title = lines[0].strip()
    num_atoms = int(lines[1].strip())
    
    atoms = []
    for i in range(2, 2 + num_atoms):
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
            'vz': 0.0,
            'raw_line': line
        }
        
        if atom['has_velocity']:
            try:
                atom['vx'] = float(line[44:52])
                atom['vy'] = float(line[52:60])
                atom['vz'] = float(line[60:68])
            except:
                atom['has_velocity'] = False
        
        atoms.append(atom)
    
    # 解析box尺寸
    box_line = lines[2 + num_atoms].strip() if len(lines) > 2 + num_atoms else ""
    
    return title, atoms, box_line


def delete_gro_atoms(atoms, atom_ids_to_delete, id_mapping):
    """
    从gro原子列表中删除指定原子，并根据id_mapping更新原子序号
    
    参数:
    atoms: 原子信息列表
    atom_ids_to_delete: 要删除的itp原子ID列表（字符串）
    id_mapping: itp旧ID到新ID的映射字典
    
    返回:
    new_atoms: 删除并重新编号后的原子列表
    """
    # 将要删除的ID转换为整数集合
    ids_to_delete = set(int(id) for id in atom_ids_to_delete)
    
    # 构建itp旧ID到新ID的映射（整数到整数）
    old_to_new = {}
    for old_id, new_id in id_mapping.items():
        old_to_new[int(old_id)] = int(new_id)
    
    new_atoms = []
    new_atom_num = 1
    
    for atom in atoms:
        old_atom_num = atom['atom_num']
        
        # 如果这个原子对应的itp原子ID在删除列表中，跳过
        if old_atom_num in ids_to_delete:
            continue
        
        # 更新原子序号
        atom['atom_num'] = old_to_new.get(old_atom_num, new_atom_num)
        new_atoms.append(atom)
    
    # 如果映射不完整，需要重新连续编号
    if len(old_to_new) < len(new_atoms):
        for i, atom in enumerate(new_atoms, 1):
            atom['atom_num'] = i
    
    return new_atoms


def renumber_gro_residues(atoms):
    """
    重新连续编号残基
    
    参数:
    atoms: 原子信息列表
    
    返回:
    atoms: 残基重新编号后的原子列表
    """
    residue_map = {}
    new_res_num = 1
    
    for atom in atoms:
        old_res_num = atom['residue_num']
        if old_res_num not in residue_map:
            residue_map[old_res_num] = new_res_num
            new_res_num += 1
        atom['residue_num'] = residue_map[old_res_num]
    
    return atoms


def write_gro_file(title, atoms, box_line, output_file):
    """
    将原子信息写入GRO文件
    
    参数:
    title: 标题行
    atoms: 原子信息列表
    box_line: box尺寸行
    output_file: 输出文件路径
    """
    with open(output_file, 'w') as f:
        # 写入标题
        f.write(f"{title}\n")
        
        # 写入原子数
        f.write(f"{len(atoms):5d}\n")
        
        # 写入原子坐标
        for atom in atoms:
            if atom['has_velocity']:
                line = (f"{atom['residue_num']:5d}"
                       f"{atom['residue_name']:5s}"
                       f"{atom['atom_name']:>5s}"
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
                       f"{atom['atom_name']:>5s}"
                       f"{atom['atom_num']:5d}"
                       f"{atom['x']:8.3f}"
                       f"{atom['y']:8.3f}"
                       f"{atom['z']:8.3f}\n")
            f.write(line)
        
        # 写入box尺寸
        f.write(f"{box_line}\n")


def process_gro_file(gro_file, atoms_to_delete, id_mapping, output_file):
    """
    处理gro文件：删除原子并重新编号
    
    参数:
    gro_file: 输入gro文件路径
    atoms_to_delete: 要删除的原子ID列表
    id_mapping: itp旧ID到新ID的映射
    output_file: 输出gro文件路径
    """
    # 解析gro文件
    title, atoms, box_line = parse_gro_file(gro_file)
    
    print(f"  原始gro文件原子数: {len(atoms)}")
    
    # 删除原子并重新编号
    new_atoms = delete_gro_atoms(atoms, atoms_to_delete, id_mapping)
    
    # 重新编号残基
    new_atoms = renumber_gro_residues(new_atoms)
    
    print(f"  删除后原子数: {len(new_atoms)}")
    print(f"  删除了 {len(atoms) - len(new_atoms)} 个原子")
    
    # 写入新文件
    write_gro_file(title, new_atoms, box_line, output_file)
    print(f"  已保存到: {output_file}")

def parse_itp_file(file_path):
    """
    解析GROMACS ITP文件, 按照[]字段分别存储内容
    """
    sections = defaultdict(list)
    current_section = None
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            
            # 跳过空行和注释行
            if not line or line.startswith(';'):
                continue
                
            # 检测新的节(section)
            section_match = re.match(r'\[(.+)\]', line)
            if section_match:
                current_section = section_match.group(1).strip()
                #print(current_section)
                continue
                
            # 将行添加到当前节
            if current_section:
                sections[current_section].append(line)
    #print(sections)
    # 处理不同节的数据
    processed_sections = {}
    
    for section_name, lines in sections.items():
        if section_name == 'moleculetype':
            processed_sections[section_name] = parse_moleculetype(lines)
        elif section_name == 'atoms':
            processed_sections[section_name] = parse_atoms(lines)
        elif section_name == 'bonds':
            processed_sections[section_name] = parse_bonds(lines)
        elif section_name == 'pairs':
            processed_sections[section_name] = parse_pairs(lines)
        elif section_name == 'angles':
            processed_sections[section_name] = parse_angles(lines)
        elif section_name == 'dihedrals':
            processed_sections[section_name] = parse_dihedrals(lines)
        else:
            processed_sections[section_name] = lines
    
    return processed_sections

def parse_moleculetype(lines):
    """解析[moleculetype]节"""
    result = {}
    if lines:
        # 第一行可能是标题行，包含';name nrexcl'
        if lines[0].startswith(';'):
            headers = [h.strip() for h in lines[0][1:].split()]
            if len(lines) > 1:
                data = lines[1].split()
                for i, header in enumerate(headers):
                    if i < len(data):
                        result[header] = data[i]
        else:
            # 没有标题行，直接解析数据
            data = lines[0].split()
            if len(data) >= 2:
                result['name'] = data[0]
                result['nrexcl'] = data[1]
    return result

def parse_atoms(lines):
    """解析[atoms]节, 按原子ID存储"""
    atoms = {}
    
    # 查找标题行
    headers = []
    for i, line in enumerate(lines):
        if line.startswith(';'):
            # 提取标题
            headers = [h.strip() for h in line[1:].split(';')[0].split()]
            break
    
    # 如果没有找到标题行，使用默认标题
    if not headers:
        headers = ['nr', 'type', 'resi', 'res', 'atom', 'cgnr', 'charge', 'mass']
    
    # 解析原子数据
    for line in lines:
        if line.startswith(';'):
            continue
            
        data = line.split()
        if len(data) < len(headers):
            continue
            
        atom_id = data[0]
        atom_info = {}
        
        for i, header in enumerate(headers):
            if i < len(data):
                # 尝试将数值转换为适当类型
                try:
                    if '.' in data[i]:
                        atom_info[header] = float(data[i])
                    else:
                        atom_info[header] = int(data[i])
                except ValueError:
                    atom_info[header] = data[i]
        
        # 添加注释信息（如果有）
        if ';' in line:
            comment = line.split(';', 1)[1].strip()
            atom_info['comment'] = comment
        
        atoms[atom_id] = atom_info
    
    return atoms

def parse_bonds(lines):
    """解析[bonds]节"""
    bonds = []
    for line in lines:
        if line.startswith(';'):
            continue
            
        data = line.split()
        if len(data) >= 3:
            bond = {
                'ai': data[0],
                'aj': data[1],
                'funct': data[2:],
                #'parameters': data[3:5],
                #'comment': ' '.join(data[5:]) if len(data) > 5 else ''
            }
            bonds.append(bond)
    
    return bonds

def parse_pairs(lines):
    """解析[pairs]节"""
    pairs = []
    for line in lines:
        if line.startswith(';'):
            continue
            
        data = line.split()
        if len(data) >= 3:
            pair = {
                'ai': data[0],
                'aj': data[1],
                'funct': data[2:],
                #'comment': ' '.join(data[3:]) if len(data) > 3 else ''
            }
            pairs.append(pair)
    
    return pairs

def parse_angles(lines):
    """解析[angles]节"""
    angles = []
    for line in lines:
        if line.startswith(';'):
            continue
            
        data = line.split()
        if len(data) >= 4:
            angle = {
                'ai': data[0],
                'aj': data[1],
                'ak': data[2],
                'funct': data[3:],
                #'parameters': data[4:6],
                #'comment': ' '.join(data[6:]) if len(data) > 6 else ''
            }
            angles.append(angle)
    
    return angles

def parse_dihedrals(lines):
    """解析[dihedrals]节"""
    dihedrals = []
    for line in lines:
        if line.startswith(';'):
            continue
            
        data = line.split()
        if len(data) >= 5:
            dihedral = {
                'ai': data[0],
                'aj': data[1],
                'ak': data[2],
                'al': data[3],
                'funct': data[4:],
                #'parameters': data[5:]
                #'phase': data[5],
                #'kd': data[6],
                #'pn': data[7],
                #'comment': ' '.join(data[8:]) if len(data) > 8 else ''
            }
            dihedrals.append(dihedral)
    
    return dihedrals

##################################################################
def delete_atoms_and_renumber(itp_data, atoms_to_delete, atom_num):
    """
    删除指定原子并重新编号所有相关部分
    
    参数:
    itp_data: 解析后的ITP数据
    atoms_to_delete: 要删除的原子ID列表
    """
    # 转换原子ID为字符串以便比较
    atoms_to_delete = [str(atom_id) for atom_id in atoms_to_delete]
    
    # 把电荷平均加到每个原子上面
    charge_delete = 0
    for delete_atom_id in atoms_to_delete:
        charge_delete += itp_data['atoms'][delete_atom_id]['charge']
    #print(charge_delete)
    avg_charge = charge_delete / atom_num

    # 1. 处理atoms节
    if 'atoms' in itp_data:
        # 删除指定原子
        for atom_id in atoms_to_delete:
            if atom_id in itp_data['atoms']:
                del itp_data['atoms'][atom_id]
        
        # 重新编号原子
        atoms = itp_data['atoms']
        sorted_atoms = sorted(atoms.items(), key=lambda x: int(x[0]))
        new_atoms = {}
        id_mapping = {}  # 存储旧ID到新ID的映射
        
        for new_id, (old_id, atom_info) in enumerate(sorted_atoms, 1):
            new_id_str = str(new_id)
            id_mapping[old_id] = new_id_str
            atom_info['nr'] = new_id
            new_atoms[new_id_str] = atom_info
            atom_info['charge'] += avg_charge

        itp_data['atoms'] = new_atoms

    # 2. 处理bonds节
    if 'bonds' in itp_data:
        new_bonds = []
        for bond in itp_data['bonds']:
            # 如果键中包含要删除的原子，则跳过
            if bond['ai'] in atoms_to_delete or bond['aj'] in atoms_to_delete:
                continue
            
            # 更新原子ID
            new_bond = bond.copy()
            new_bond['ai'] = id_mapping.get(bond['ai'], bond['ai'])
            new_bond['aj'] = id_mapping.get(bond['aj'], bond['aj'])
            new_bonds.append(new_bond)
        
        itp_data['bonds'] = new_bonds
    
    # 3. 处理pairs节
    if 'pairs' in itp_data:
        new_pairs = []
        for pair in itp_data['pairs']:
            # 如果对中包含要删除的原子，则跳过
            if pair['ai'] in atoms_to_delete or pair['aj'] in atoms_to_delete:
                continue
            
            # 更新原子ID
            new_pair = pair.copy()
            new_pair['ai'] = id_mapping.get(pair['ai'], pair['ai'])
            new_pair['aj'] = id_mapping.get(pair['aj'], pair['aj'])
            new_pairs.append(new_pair)
        
        itp_data['pairs'] = new_pairs
    
    # 4. 处理angles节
    if 'angles' in itp_data:
        new_angles = []
        for angle in itp_data['angles']:
            # 如果角中包含要删除的原子，则跳过
            if (angle['ai'] in atoms_to_delete or 
                angle['aj'] in atoms_to_delete or 
                angle['ak'] in atoms_to_delete):
                continue
            
            # 更新原子ID
            new_angle = angle.copy()
            new_angle['ai'] = id_mapping.get(angle['ai'], angle['ai'])
            new_angle['aj'] = id_mapping.get(angle['aj'], angle['aj'])
            new_angle['ak'] = id_mapping.get(angle['ak'], angle['ak'])
            new_angles.append(new_angle)
        
        itp_data['angles'] = new_angles
    
    # 5. 处理dihedrals节
    if 'dihedrals' in itp_data:
        new_dihedrals = []
        for dihedral in itp_data['dihedrals']:
            # 如果二面角中包含要删除的原子，则跳过
            if (dihedral['ai'] in atoms_to_delete or 
                dihedral['aj'] in atoms_to_delete or 
                dihedral['ak'] in atoms_to_delete or 
                dihedral['al'] in atoms_to_delete):
                continue
            
            # 更新原子ID
            new_dihedral = dihedral.copy()
            new_dihedral['ai'] = id_mapping.get(dihedral['ai'], dihedral['ai'])
            new_dihedral['aj'] = id_mapping.get(dihedral['aj'], dihedral['aj'])
            new_dihedral['ak'] = id_mapping.get(dihedral['ak'], dihedral['ak'])
            new_dihedral['al'] = id_mapping.get(dihedral['al'], dihedral['al'])
            new_dihedrals.append(new_dihedral)
        
        itp_data['dihedrals'] = new_dihedrals
    
    return itp_data, id_mapping

def write_itp_file(itp_data, output_file):
    """
    将ITP数据写入文件
    """
    with open(output_file, 'w+') as f:
        # 写入moleculetype节
        if 'moleculetype' in itp_data:
            f.write('[ moleculetype ]\n')
            moltype = itp_data['moleculetype']
            f.write(';name            nrexcl\n')
            #print(moltype['name'])
            f.write("%s     %s\n" % (moltype['name'], moltype['nrexcl']))
            f.write('\n')

        # 写入atoms节
        if 'atoms' in itp_data:
            f.write('[ atoms ]\n')
            f.write(';   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type\n')
            
            atoms = itp_data['atoms']
            for atom_id in sorted(atoms.keys(), key=int):
                atom = atoms[atom_id]
                line = f"{atom['nr']:>5} {atom['type']:>5} {atom['resi']:>5} {atom['res']:>4} {atom['atom']:>5} {atom['cgnr']:>5} {atom['charge']:>10.6f} {atom['mass']:>10.5f}"
                if 'comment' in atom:
                    line += f" ; {atom['comment']}"
                f.write(line + '\n')
            f.write('\n')
        
        # 写入bonds节
        if 'bonds' in itp_data and itp_data['bonds']:
            f.write('[ bonds ]\n')
            f.write(';   ai     aj funct   r             k\n')
            for bond in itp_data['bonds']:
                line = f"{bond['ai']:>6} {bond['aj']:>6} {' '.join(str(p) for p in bond['funct']):>7}"
                if bond.get('comment'):
                    line += f" ; {bond['comment']}"
                f.write(line + '\n')
            f.write('\n')
        
        # 写入pairs节
        if 'pairs' in itp_data and itp_data['pairs']:
            f.write('[ pairs ]\n')
            f.write(';   ai     aj    funct\n')
            for pair in itp_data['pairs']:
                line = f"{pair['ai']:>6} {pair['aj']:>6} {' '.join(str(p) for p in pair['funct']):>7}"
                if pair.get('comment'):
                    line += f" {pair['comment']}"
                f.write(line + '\n')
            f.write('\n')
        
        # 写入angles节
        if 'angles' in itp_data and itp_data['angles']:
            f.write('[ angles ]\n')
            f.write(';   ai     aj     ak    funct   theta         cth\n')
            for angle in itp_data['angles']:
                line = f"{angle['ai']:>6} {angle['aj']:>6} {angle['ak']:>6} {' '.join(str(p) for p in angle['funct']):>7}"
                if angle.get('comment'):
                    line += f" {angle['comment']}"
                f.write(line + '\n')
            f.write('\n')
        
        # 写入dihedrals节
        if 'dihedrals' in itp_data and itp_data['dihedrals']:
            f.write('[ dihedrals ]\n')
            f.write(';   ai     aj     ak     al   func   phase     kd      pn\n')
            for dihedral in itp_data['dihedrals']:
                line = f"{dihedral['ai']:>6} {dihedral['aj']:>6} {dihedral['ak']:>6} {dihedral['al']:>6}  {' '.join(str(p) for p in dihedral['funct']):>7}"
                if dihedral.get('comment'):
                    line += f" {dihedral['comment']}"
                f.write(line + '\n')
            f.write('\n')
        '''
        # 写入其他节
        for section_name in itp_data:
            if section_name not in ['moleculetype', 'atoms', 'bonds', 'pairs', 'angles', 'dihedrals']:
                f.write(f'[ {section_name} ]\n')
                for line in itp_data[section_name]:
                    f.write(line + '\n')
                f.write('\n')
        '''
# 使用示例
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='GROMACS ITP/GRO 原子删除工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 只处理ITP文件
  python delete_itp_atom.py -i TMC.itp -d 11 13 15
  
  # 同时处理ITP和GRO文件
  python delete_itp_atom.py -i TMC.itp -g TMC.gro -d 11 13 15
  
  # 指定输出文件名
  python delete_itp_atom.py -i TMC.itp -g TMC.gro -d 11 13 15 -oi TMC_new.itp -og TMC_new.gro
        """
    )
    
    parser.add_argument('-i', '--itp', required=True,
                        help='输入的ITP文件路径')
    parser.add_argument('-g', '--gro', default=None,
                        help='输入的GRO文件路径（可选）')
    parser.add_argument('-d', '--delete', nargs='+', type=int, required=True,
                        help='要删除的原子ID列表')
    parser.add_argument('-oi', '--output-itp', default=None,
                        help='输出的ITP文件路径（默认: 输入文件名_modified.itp）')
    parser.add_argument('-og', '--output-gro', default=None,
                        help='输出的GRO文件路径（默认: 输入文件名_modified.gro）')
    
    args = parser.parse_args()
    
    # 设置默认输出文件名
    itp_base = os.path.splitext(args.itp)[0]
    output_itp = args.output_itp if args.output_itp else f"{itp_base}_modified.itp"
    output_gro = args.output_gro if args.output_gro else (f"{os.path.splitext(args.gro)[0]}_modified.gro" if args.gro else None)
    
    itp_file = args.itp
    gro_file = args.gro
    atoms_to_delete = args.delete
    
    print("=" * 60)
    print("GROMACS ITP/GRO 原子删除工具")
    print("=" * 60)
    
    # 检查输入文件是否存在
    if not os.path.exists(itp_file):
        print(f"错误: ITP文件不存在: {itp_file}")
        exit(1)
    
    if gro_file and not os.path.exists(gro_file):
        print(f"错误: GRO文件不存在: {gro_file}")
        exit(1)
    
    # 处理ITP文件
    print(f"\n[1] 处理ITP文件: {itp_file}")
    print(f"    要删除的原子ID: {atoms_to_delete}")
    
    itp_data = parse_itp_file(itp_file)
    atoms = itp_data.get('atoms')
    print(f"    原始原子数: {len(atoms)}")
    
    atom_num = len(atoms) - len(atoms_to_delete)
    print(f"    删除后原子数: {atom_num}")
    
    # 打印要删除的原子信息
    print(f"\n    要删除的原子:")
    for atom_id in atoms_to_delete:
        atom_id_str = str(atom_id)
        if atom_id_str in atoms:
            atom_info = atoms[atom_id_str]
            print(f"      ID={atom_id}: {atom_info['atom']} ({atom_info['type']})")
        else:
            print(f"      ID={atom_id}: 未找到")
    
    # 删除原子并重新编号
    modified_data, id_mapping = delete_atoms_and_renumber(itp_data, atoms_to_delete, atom_num)
    write_itp_file(modified_data, output_itp)
    print(f"\n    已保存修改后的ITP文件: {output_itp}")
    
    # 打印ID映射
    print(f"\n    原子ID映射 (旧 -> 新):")
    for old_id, new_id in sorted(id_mapping.items(), key=lambda x: int(x[0])):
        print(f"      {old_id} -> {new_id}")
    
    # 处理GRO文件
    if gro_file:
        print(f"\n[2] 处理GRO文件: {gro_file}")
        process_gro_file(gro_file, atoms_to_delete, id_mapping, output_gro)
    
    print("\n" + "=" * 60)
    print("处理完成!")
    print("=" * 60)