# -*- coding: utf-8 -*-
import re
from collections import defaultdict

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
                'funct': data[2],
                'parameters': data[3:5],
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
                'funct': data[2],
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
        if len(data) >= 3:
            angle = {
                'ai': data[0],
                'aj': data[1],
                'ak': data[2],
                'funct': data[3],
                'parameters': data[4:6],
                'comment': ' '.join(data[6:]) if len(data) > 6 else ''
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
                'funct': data[4],
                'parameters': data[5:]
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
                if len(bond['parameters']) == 0:
                    line = f"{bond['ai']:>6} {bond['aj']:>6} {bond['funct']:>5}"
                else:
                    line = f"{bond['ai']:>6} {bond['aj']:>6} {bond['funct']:>5} {bond['parameters'][0]:>12} {bond['parameters'][1]:>12}"
                if bond.get('comment'):
                    line += f" ; {bond['comment']}"
                f.write(line + '\n')
            f.write('\n')
        
        # 写入pairs节
        if 'pairs' in itp_data and itp_data['pairs']:
            f.write('[ pairs ]\n')
            f.write(';   ai     aj    funct\n')
            for pair in itp_data['pairs']:
                line = f"{pair['ai']:>6} {pair['aj']:>6} {pair['funct']:>7}"
                if pair.get('comment'):
                    line += f" {pair['comment']}"
                f.write(line + '\n')
            f.write('\n')
        
        # 写入angles节
        if 'angles' in itp_data and itp_data['angles']:
            f.write('[ angles ]\n')
            f.write(';   ai     aj     ak    funct   theta         cth\n')
            for angle in itp_data['angles']:
                if len(angle['parameters']) == 0:
                    line = f"{angle['ai']:>6} {angle['aj']:>6} {angle['ak']:>6} {angle['funct']:>7}"
                else:
                    line = f"{angle['ai']:>6} {angle['aj']:>6} {angle['ak']:>6} {angle['funct']:>7} {angle['parameters'][0]:>12} {angle['parameters'][1]:>12}"
                if angle.get('comment'):
                    line += f" {angle['comment']}"
                f.write(line + '\n')
            f.write('\n')
        
        # 写入dihedrals节
        if 'dihedrals' in itp_data and itp_data['dihedrals']:
            f.write('[ dihedrals ]\n')
            f.write(';   ai     aj     ak     al   func   phase     kd      pn\n')
            for dihedral in itp_data['dihedrals']:
                if len(dihedral['parameters']) == 0:
                    line = f"{dihedral['ai']:>6} {dihedral['aj']:>6} {dihedral['ak']:>6} {dihedral['al']:>6} {dihedral['funct']:>5}"
                else:
                    line = f"{dihedral['ai']:>6} {dihedral['aj']:>6} {dihedral['ak']:>6} {dihedral['al']:>6} {dihedral['funct']:>5} {' '.join(str(p) for p in dihedral['parameters']):>7}"
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
    #######################################################
    # itp 名称
    itp_file = "protein.itp"
    # 要删除的原子ID列表
    atoms_to_delete = [15]  # 

    #######################################################
    itp_data = parse_itp_file(itp_file)
    #print(itp_data)
    moleculetype = itp_data.get('moleculetype')
    #print(moleculetype)
    atoms = itp_data.get('atoms')
    #print(atoms)
    bonds = itp_data.get('bonds')
    #print(bonds)
    pairs = itp_data.get('pairs')
    #print(pairs)
    angles = itp_data.get('angles')
    #print(angles)
    dihedrals = itp_data.get('dihedrals')
    #print(dihedrals)
    atom_num = len(atoms) - len(atoms_to_delete)
    #print(atom_num)
    modified_data, id_mapping = delete_atoms_and_renumber(itp_data, atoms_to_delete, atom_num)
    write_itp_file(modified_data, "modified.itp")