#!/usr/bin/env python3

import linecache
import os
import sys
import math
import numpy as np
import re
import time
import datetime
import subprocess
from collections import defaultdict
from itertools import combinations
from itertools import product
import itertools

def ReadGMXGro(GroFile):
    '''
    Read single frame .gro file,
    Return box, moles. moles: list consisted of Molecule class.
    '''
    with open(GroFile, 'r') as f:
        lines = [line.rstrip() for line in f.readlines()]
    totalmoles = []
    openedfile = open(GroFile,'r')
    natoms = len(openedfile.readlines()) - 3
    boxv = np.array(list(map(float, lines[2 + natoms].split())))

    atoms = lines[2: 3 + natoms]
    n = 0
    while n < natoms:
        name = atoms[n][:10].strip()
        typename = atoms[n][5:10].strip()
        coordinates = []
        symbols = []
        mole = []
        for nextn in range(n, natoms + 1):
            coordinates.append(list(map(float, atoms[nextn][20:44].split())))
            symbols.append((atoms[nextn][10:15]).strip())
            if nextn + 1 == natoms or atoms[nextn + 1][:10].strip() != name:
                break
        line1 = 0
        for info1 in symbols:
            atomlist = []
            x = coordinates[line1][0]
            y = coordinates[line1][1]
            z = coordinates[line1][2]
            atomlist.append(typename)
            atomlist.append(info1)
            atomlist.append(x)
            atomlist.append(y)
            atomlist.append(z)
            mole.append(atomlist)
            line1 += 1
        totalmoles.append(mole)
        n = nextn + 1
    openedfile.close()
    f.close()
    return natoms, boxv, totalmoles
#
def WirteGMXGro(totalmoles, NewGroName, natoms, boxsize):
    '''
    Read totalmoles: list
    Rcreate new gmx gro file.
    '''
    NewFile = open(NewGroName, 'w+')
    NewFile.write("New Gro File\n")
    NewFile.write("%5d\n" % natoms)
    n = 1
    m = 1
    #read total molecle
    if len(totalmoles) > 0:
        for info1 in totalmoles:
            #print(info1)
            resid = n
            # read single molecle
            for info2 in info1:
                # read single atomic information
                resname = info2[0]
                atomtype = info2[1]
                x = info2[2]
                y = info2[3]
                z = info2[4]
                natom = (m - 1) % 99999 + 1
                NewFile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %
                              (resid, resname, atomtype, natom, x, y, z))
                m += 1
            n += 1
    else:
        resid = 0
        natom = 0
    NewFile.write("%10.5f%10.5f%10.5f\n" % (boxsize[0], boxsize[1], boxsize[2]))
    NewFile.flush()
    NewFile.close()
    return resid, natom

def GetCentCoord(atom_list):
    x_sum = 0
    y_sum = 0
    z_sum = 0
    com_coord = []
    for info in atom_list:
        x = info[-3]
        y = info[-2]
        z = info[-1]
        x_sum += float(x)
        y_sum += float(y)
        z_sum += float(z)

    x_com = round(x_sum / len(atom_list),3)
    y_com = round(y_sum / len(atom_list),3)
    z_com = round(z_sum / len(atom_list),3)

    com_coord.append(x_com)
    com_coord.append(y_com)
    com_coord.append(z_com)

    return com_coord

def GetItpFrag(FileName, FragNameBegin, FragNameEnd="",Order=0):
    '''
    get itp file fragement infromation, [ bonds] or [ angle ] etc.
    :param FileName: itp file name
    :param FragNameBegin,FragNameEnd: bonds  angle  pair  etc.
    :return: FragmentList list
    '''
    FragmentList = []
    TwoDimensionList = []
    keyword1 = r'\s?\[\s.*%s\s.*\]\s.*' % FragNameBegin
    keyword2 = r'\s?\[\s.*%s\s.*\]\s.*' % FragNameEnd
    LineBegin, TotalLine = GetLineNum(FileName, keyword1)
    LineBegin = LineBegin[Order] + 1
    while LineBegin <= TotalLine:
        LineData = str(linecache.getline(FileName, LineBegin))
        if re.match(keyword2, LineData) != None:
            break
        elif re.match(r'\s?;\s?.*', LineData) == None and len(LineData.strip()) > 1:
            FragmentList.append(LineData)
            TwoDimensionList.append(LineData.split())
        LineBegin += 1
    linecache.clearcache()
    return FragmentList,TwoDimensionList

def GetLineNum(FileName, KeyWord):
    '''
    fine line number
    :param FileName:
    :param KeyWord:
    :return: LineList list and LineTotal
    '''
    filename = open(FileName,'r')
    LineTotal = len(filename.readlines())
    n = 0
    LineList = []
    while n <= LineTotal:
        LineData = str(linecache.getline(FileName, n))
        #info = LineData.find(KeyWord)
        #KeyWord1 =  re.compile(r'\s?\[\s.*%s\s.*\]\s.*' % KeyWord)
        #print(KeyWord)
        info = re.search(KeyWord, LineData)
        #info = KeyWord1.search(LineData)
        if info == None:
            pass
        else:
            LineList.append(n)
        n += 1
    filename.flush()
    filename.close()
    return LineList, LineTotal
#
def GetFile(FileName):
    '''
    Read file as a list
    :param FileName:
    :return: FileDetaill list
    '''
    FileDetail = []
    with open(FileName) as file:
        for info in file:
            FileDetail.append(info)
    file.flush()
    file.close()
    return FileDetail

def AllItpfragDict(FileName,TopParentList):
    '''
    split file to fragment of  ParentList
    :param FileName:
    :param TopParentList:
    :return: Dict of fragment about all parentlist
    '''
    FragDict = {}
    for info in TopParentList:
        LineNum = GetLineNum(FileName, info)[0]
        LineLen = len(LineNum)
        if LineLen == 1:
            FragmentList,TwoDimensionList = GetItpFrag(FileName, info, "",LineLen-1)
            FragDict.update({info: TwoDimensionList})
        else:
            m = 1
            TmpList = []
            while m <= LineLen:
                FragmentList,TwoDimensionList = GetItpFrag(FileName, info, "",m-1)
                TmpList += TwoDimensionList
                m += 1
            FragDict.update({info: TmpList})
    return FragDict

def WriteNewItp(new_itp_file, Keywrod, fragment, res_name, new_atoms, terminal_atoms, origin_atom_id, AddNum=0):
    if Keywrod == "atoms":
        if  origin_atom_id <= len(fragment):
            #print(len(fragment))
            #print(fragment[origin_atom_id - 1][4])
            old_name = fragment[origin_atom_id - 1][4]
            if old_name in terminal_atoms:
                index = terminal_atoms.index(old_name)
                new_atom = new_atoms[index]
                #print(new_atom)
                fragment[origin_atom_id - 1][4] = new_atom
        for info in fragment:
            if info[3] == res_name:
                new_resname = info[3]
            else:
                new_resname = res_name
            new_itp_file.write("%6d%5s%6d%6s%6s%5d%13.5f%13.5f\n" %
                            (int(info[0]) + AddNum, info[1], int(info[2]), new_resname, info[4],
                             int(info[5]) + AddNum, float(info[6]), float(info[7])))
    elif Keywrod == "bonds":
        for info in fragment:
            new_itp_file.write("%6d%7d%4d" % (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])))
            n = 1
            for info1 in info:
                if n >= 4:
                    new_itp_file.write("   %s" % info1)
                n += 1
            new_itp_file.write("\n")
    elif Keywrod == "pairs":
        for info in fragment:
            new_itp_file.write("%6d%7d%4d" % (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])))
            n = 1
            for info1 in info:
                if n >= 4:
                    new_itp_file.write("   %s" % info1)
                n += 1
            new_itp_file.write("\n")
    elif Keywrod == "angles":
        for info in fragment:
            new_itp_file.write("%6d%7d%7d%7d" % (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,int(info[3])))
            n = 1
            for info1 in info:
                if n >= 5:
                    new_itp_file.write("   %s" % info1)
                n += 1
            new_itp_file.write("\n")
    elif Keywrod == "dihedrals":
        for info in fragment:
            new_itp_file.write("%6d%7d%7d%7d%7d" % (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                                     int(info[3])+AddNum,int(info[4])))
            n = 1
            for info1 in info:
                if n >= 6:
                    new_itp_file.write("   %s" % info1)
                n += 1
            new_itp_file.write("\n")
    else:
        print("!!! Unrecognized keyword: %s !!!" % Keywrod)
        sys.exit()

def CombineITP(mol_atom_num, repeat_num, new_itp_name, monomer_itp_name, all_bond_info, res_name, new_atoms, terminal_atoms, origin_atom_id):
    TopParentList = ["atoms", "bonds", "pairs", "angles", "dihedrals"]
        # 
    args = []
    for info in repeat_num:
        args.append(monomer_itp_name[info])
        args.append('1')
    if (len(args) % 2) == 0:
        pass
    else:
        print("Incorrect input,Please Check!!!")
        sys.exit()
    itp_name = []
    itp_name_num = {}
    n1 = 1
    for info in args:
        if (n1 % 2) != 0:
            itp_name.append(args[n1-1])
            itp_name_num[args[n1-1]] = args[n1]
        else:
            pass
        n1 += 1
    #print(itp_name)
    new_itp_file = open(new_itp_name+".itp", 'w+')
    new_itp_file.write("; Merge itp files (%s)\n" % ', '.join(itp_name))
    new_itp_file.write("[ moleculetype ]\n; Name   nrexcl\n%s     3\n" % res_name)
    all_itp_info = []
    all_atom_mol_num = []
    for info in itp_name:
        atom_num = len(GetItpFrag(info, "atoms", FragNameEnd="", Order=0)[0])
        mol_num = int(itp_name_num[info])
        frag_dict = AllItpfragDict(info,TopParentList)
        all_itp_info.append(frag_dict)
        all_atom_mol_num.append(atom_num)
        all_atom_mol_num.append(mol_num)
    for info1 in TopParentList:
        new_itp_file.write("[ %s ]\n" % info1)
        n2 = 1
        #n4 = -int(all_atom_mol_num[0])
        n4 = 0
        for info2 in all_itp_info:
            n3 = 1
            mol_num1 = int(all_atom_mol_num[n2])
            atom_num1 = int(all_atom_mol_num[n2 - 1])
            while n3 <= mol_num1:
                #n4 += atom_num1
                #print(mol_num1,atom_num1)
                #print(n4)
                fragment = info2[info1]
                WriteNewItp(new_itp_file, info1, fragment, res_name, new_atoms, terminal_atoms, origin_atom_id, n4)
                n4 += atom_num1
                n3 += 1
            n2 += 2
    # 添加链接关系
    total_mol_bonded_info = all_bond_info
    #print(total_mol_bonded_info["10"]["bonds"])
    mol_num = len(repeat_num)
    atom1_num = 0
    atom2_num = 0
    for mol_id in range(mol_num-1):
        id_b = repeat_num[mol_id]
        id_e = repeat_num[mol_id+1]
        if str(id_b) == '0' and str(id_e) == '1':
            id_str = str(id_b)+str(id_e)
        else:
            continue
        mol_bonded_info = total_mol_bonded_info[id_str]
        atom2_num += mol_atom_num[id_b]
        #print(mol_bonded_info)
        for info in mol_bonded_info["bonds"]:
            #print(info)
            atom1 = int(info[0])
            atom2 = int(info[1])
            funct = int(info[2])
            r0 = info[3]
            k = info[4]
            #print(atom1,atom2)
            #other = ""
            #for oth in info[2:]:
            #    other += "%15.10f" % oth
            ##print(other)
            atom1 += atom1_num
            atom2 += atom2_num 
            #print(atom1, atom2)
            #atom1_num = atom2_num
            bonds_info = " %6d" % atom1 + " %7d" % atom2 + " %4d" % funct + " %15.6f" % r0 + " %15.6f" % k + " ; user define"
            new_itp_file.write(" [ bonds ]\n")
            new_itp_file.write("%s\n" % bonds_info)
            
        for info in mol_bonded_info["angles1"]:
            atom1 = int(info[0])
            atom2 = int(info[1])
            atom3 = int(info[2])
            funct = int(info[3])
            theta0 = float(info[4])
            k = float(info[5])
            #print(other)
            atom1 += atom1_num
            atom2 += atom1_num
            atom3 += atom2_num
            angles_info = " %6d" % atom1 + " %7d" % atom2 + " %7d" % atom3 + " %4d" % funct + " %15.6f" % theta0 + " %15.6f" % k + " ; user define"
            #print(angles_info)
            new_itp_file.write(" [ angles ]\n")
            new_itp_file.write("%s\n" % angles_info)
        for info in mol_bonded_info["angles2"]:
            atom1 = int(info[0])
            atom2 = int(info[1])
            atom3 = int(info[2])
            funct = int(info[3])
            theta0 = float(info[4])
            k = float(info[5])
            #print(other)
            atom1 += atom1_num
            atom2 += atom2_num
            atom3 += atom2_num
            angles_info = " %6d" % atom1 + " %7d" % atom2 + " %7d" % atom3 + " %4d" % funct + " %15.6f" % theta0 + " %15.6f" % k + " ; user define"
            new_itp_file.write(" [ angles ]\n")
            new_itp_file.write("%s\n" % angles_info)  
        for info in mol_bonded_info["dihedrals"]:
            atom1 = int(info[0])
            atom2 = int(info[1])
            atom3 = int(info[2])
            atom4 = int(info[3])
            other = ""
            for oth in info[4:]:
                other += "%15s" % oth
            #print(other)
            atom1 += atom1_num
            atom2 += atom1_num
            atom3 += atom2_num
            atom4 += atom2_num
            dihedral_info = " %10d" % atom1 + " %10s" % atom2 + " %10s" % atom3 + " %10s" % atom4 + other + " ; user define"
            #print(dihedral_info)
            new_itp_file.write(" [ dihedrals ]\n")
            new_itp_file.write("%s\n" % dihedral_info) 
            
        atom1_num = atom2_num
        
    new_itp_file.flush()
    new_itp_file.close()
    return

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
##################################################################
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
##################################################################
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
##################################################################
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
        file.flush()
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

def ReadBonded(file_name):
    file_list = []
    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()  
            if line == '' or ';' in line:
                continue
            else:
                words = line.split()
                file_list.append(words)
        file.flush()
    bonded_dict = {}
    bonded_list = []
    num = 0
    for info in file_list:
        num += 1
        if info[0] == "[":
            if len(bonded_list) > 0:
                bonded_dict[bonded_info] = bonded_list
            bonded_info = info[1]
            bonded_list = []
        if info[0] != "[":
            bonded_list.append(info)
        if num == len(file_list):
            bonded_dict[bonded_info] = bonded_list
    return bonded_dict

def parse_itp_bonds(itp_file, target_atom_id):
    """
    解析ITP文件，找到与指定原子有键连接的所有原子
    
    参数:
        itp_file: ITP文件路径
        target_atom_id: 目标原子ID
    
    返回:
        二维列表: [[id1, type1], [id2, type2], ...]
    """
    atoms = {}  # 原子ID -> 类型
    bonds = []  # 键对列表
    
    section = None
    
    with open(itp_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # 检测section开始
            if line.startswith('['):
                match = re.match(r'\[(\s*)?(\w+)(\s*)?\]', line)
                if match:
                    section = match.group(2).lower()
                continue
            
            # 跳过空行和注释
            if not line or line.startswith(';'):
                continue
            
            # 解析atoms部分
            if section == 'atoms':
                # 使用正则表达式匹配数字开头
                if re.match(r'^\d', line):
                    parts = line.split()
                    try:
                        atom_id = int(parts[0])
                        atom_type = parts[1]
                        atoms[atom_id] = atom_type
                    except (ValueError, IndexError):
                        continue
            
            # 解析bonds部分
            elif section == 'bonds':
                # 使用正则表达式匹配数字开头
                if re.match(r'^\d', line):
                    parts = line.split()
                    try:
                        atom1 = int(parts[0])
                        atom2 = int(parts[1])
                        bonds.append((atom1, atom2))
                    except (ValueError, IndexError):
                        continue
        f.flush()
    # 查找键连接的原子
    result = []
    for a1, a2 in bonds:
        if a1 == target_atom_id and a2 in atoms:
            result.append([a2, atoms[a2]])
        elif a2 == target_atom_id and a1 in atoms:
            result.append([a1, atoms[a1]])
    
    return result, [target_atom_id, atoms[target_atom_id]]

def RunGmx(TopFile,Number,MinimFile,MdFile):
    '''
    Run Gmx program
    '''
    print("Gromacs program is runing...")
    # call gromacs program
    os.system('gmx grompp -f %s -c md%s.gro -p %s -o em%s.tpr -maxwarn 12 >& pre_em.log' % (MinimFile,Number,TopFile,Number))
    os.system('gmx mdrun -deffnm em%s -nt 6 -v >& run_em.log' % Number)
    os.system('gmx grompp -f %s -c em%s.gro -p %s -o md%s.tpr -maxwarn 12 >& pre_md.log' % (MdFile,Number,TopFile,Number))
    os.system('gmx mdrun -deffnm md%s  -nt 6 -v >& run_md.log' % Number)
    #os.system('echo 0 | gmx trjconv -f md%s.gro -s md%s.tpr -o md%s.gro -pbc whole >& /dev/null' % (Number,Number,Number))
    os.system('rm  ./#*  *.edr *.trr *.tpr *.cpt *.xtc em*.gro md*.log em*.log >& /dev/null')

def main(em_mdp, md_mdp, top_file, gro_file, terminal_atoms, new_atoms, link_atoms, terminal_gros, terminal_itps, bonded_info, iter_bond):

    # 把封端的gro内容全部写入 all_terminal_gros 字典
    all_terminal_gros = {}
    all_termianl_CentCoord = {}
    for index, value in enumerate(terminal_atoms):
        totalmoles = ReadGMXGro(terminal_gros[index])[2]
        all_terminal_gros[value] = totalmoles[0]
        all_termianl_CentCoord[value] = GetCentCoord(totalmoles[0])
    #print(all_terminal_gros["A1"])
    #print(all_termianl_CentCoord["A1"])
    #print(all_termianl_CentCoord["D1"])
    #print(totalmoles)

    result = subprocess.run(
        "grep -o ' %s ' init.gro | wc -l" % terminal_atoms[0],
        shell=True,                # 启用 shell 以便使用管道
        capture_output=True,       # 捕获标准输出和标准错误
        text=True                  # 以文本形式返回（而非字节）
    )
    atom1_num = int(result.stdout.strip())
    #print(atom1_num)
    result = subprocess.run(
        "grep -o ' %s ' %s | wc -l" % (terminal_atoms[1], gro_file),
        shell=True,                # 启用 shell 以便使用管道
        capture_output=True,       # 捕获标准输出和标准错误
        text=True                  # 以文本形式返回（而非字节）
    )
    atom2_num = int(result.stdout.strip())
    terminal_atom_total = atom1_num + atom2_num
    #print(terminal_atom_total)
    total_iter = math.ceil(terminal_atom_total/iter_bond)
    iter = 1
    while iter <= total_iter:  # 循环添加封端
        
        if iter == 1:
            gro = gro_file
        else:
            gro = "md%s.gro" % (iter - 1)
        #print(gro)
        # 读取gro文件
        natoms = ReadGMXGro(gro)[0]
        boxsize = ReadGMXGro(gro)[1]
        totalmoles = ReadGMXGro(gro)[2]
        # 标记出具有A1和D1的分子顺序
        A1_D1_mol_id = []
        for mol_index, mol_info in enumerate(totalmoles):
            for atom_index, atom_info in enumerate(mol_info):
                if atom_info[1] in terminal_atoms:
                    A1_D1_mol_id.append([mol_index, atom_index, atom_info[1]])
        #print(A1_D1_mol_id)
        tmp_n = 0
        add_atom_num = 0
        if len(A1_D1_mol_id) < iter_bond:
            iter_bond = len(A1_D1_mol_id)
        while tmp_n < iter_bond: # 每次添加iter_bond数量封端
            #print(tmp_n)
            # 识别封端的信息
            mol_id = A1_D1_mol_id[tmp_n][0]  # gro的第几个分子
            #print("1111")
            #print(mol_id)
            atom_id = A1_D1_mol_id[tmp_n][1] # 分子中第几个原子
            #print(atom_id)
            #sys.exit()
            terminal_atom_info = totalmoles[mol_id][atom_id]
            res_name = terminal_atom_info[0]
            atom_name = terminal_atom_info[1]
            # 识别移动封端gro
            move_mole = all_terminal_gros[atom_name]
            # 修改残基名
            for i, item in enumerate(move_mole):
                move_mole[i][0] = res_name
            # 识别到A1之后,把OH.gro移动到A1旁边,并生成totalmoles形式的列表,下面以不移动的代替
            #print(all_termianl_CentCoord["A1"])
            # 计算A1/D1到封端原子的向量
            #print(totalmoles[mol_id][atom_id])
            vector = np.array(totalmoles[mol_id][atom_id][2:]) - np.array(all_termianl_CentCoord[atom_name]) - 0.1  # 靠近目标A1/D1距离0.1nm
            after_moved = []
            for move_info in move_mole:
                #print(move_info)
                res_move = move_info[0]
                atom_move = move_info[1]
                x_move = move_info[2] + vector[0]
                y_move = move_info[3] + vector[1]
                z_move = move_info[4] + vector[2]
                after_moved.append([res_move, atom_move, x_move, y_move, z_move])
            #print(move_mole)
            #print(after_moved)
            #sys.exit()
            # 把移动后的gro添加到总gro里面
            #print(totalmoles[mol_id])
            for add_info in after_moved:
                totalmoles[mol_id].append(add_info)
                # 替换封端后A1/D1原子的名称
                atom_index = terminal_atoms.index(atom_name)
                totalmoles[mol_id][atom_id][1] = new_atoms[atom_index]
                #print(totalmoles[mol_id][atom_id][1])
                #print(add_info)
                add_atom_num += 1
            #print(totalmoles[mol_id])
            # 针对res_name的itp进行添加
            origin_itp = "%s.itp" %  res_name
            # 被添加的itp里面链接原子的id
            origin_atom_id = atom_id + 1
            #print(origin_atom_id)
            # 添加的itp名称和原子序号
            add_itp = terminal_itps[atom_index]
            # 添加的itp里面链接原子id
            add_atom_id = link_atoms[atom_index]
            #print(add_atom_id)
            #print(add_itp)
            #print(origin_itp)
            # 1. 添加合并后的itp的成键关系，读取gaff_bonded.dat
            monomer_itp_name = [origin_itp, add_itp]
            link_atom_id_modified = [[origin_atom_id, origin_atom_id], [add_atom_id, add_atom_id]]
            repeat_num = [0, 1]
            origin_itp_info = parse_itp_file(origin_itp)
            add_itp_info = parse_itp_file(add_itp)
            mol_atom_num = [len(origin_itp_info["atoms"]), len(add_itp_info["atoms"])]
            #print(mol_atom_num)
            #print(link_atom_id_modified)
            # 识别和链接原子成键的原子的id和type，以及链接原子自身的type
            bond_info = []  # 和链接原子成键的原子的id和type
            link_info = []  # 链接原子自身的type
            for index, value in enumerate(link_atom_id_modified):
                tmp1 = []
                tmp2 = []
                for link in value:
                    #print(monomer_itp_name[index])
                    link1, link2 = parse_itp_bonds("%s" % monomer_itp_name[index], link)
                    #print(link1, link2)
                    tmp1.append(link1)
                    tmp2.append(link2)
                bond_info.append(tmp1)
                link_info.append(tmp2)
            #print(bond_info)
            #print(link_info)
            bond_param = ReadBonded("gaff_bonded.dat")
            #print(file['bonds'])
            all_bond_info = {}
            mol_list = ["%s" % i for i in range(len(mol_atom_num))]
            #print(mol_list)
            for mol_pair in itertools.product(mol_list, repeat=2):
                # 循环每个分子组合顺序
                pair_bond_info = {}
                if mol_pair[0] == '0' and mol_pair[1] == '1':
                    pair = mol_pair[0]+mol_pair[1]
                else:
                    continue
                mol_id_1 = int(mol_pair[0])
                mol_id_2 = int(mol_pair[1])
                #print(mol_id_1)
                # 添加bond
                atom_head_id = str(link_info[mol_id_1][1][0])
                atom_tail_id = str(link_info[mol_id_2][0][0])
                atom_head_type = link_info[mol_id_1][1][1]
                atom_tail_type = link_info[mol_id_2][0][1]
                # 插入寻找type链接关系的函数，返回r0和k
                bonds_p = bond_param['bonds']
                found_param = False
                for param in bonds_p:
                    bond_type1 = atom_head_type + '-' + atom_tail_type
                    bond_type2 = atom_tail_type + '-' + atom_head_type
                    if param[0] == bond_type1 or param[0] == bond_type2:
                        k = float(param[1]) * 4.184 * 100 * 2  # kcal/mol/A^2 to kJ/mol/nm^2
                        r0 = float(param[2]) * 0.1 #  A to nm
                        found_param = True
                if found_param == False:
                    print("!!!! No parameter %s or %s was found. Please check the bonged.dat file." % (bond_type1, bond_type2))
                    sys.exit()
                pair_bond_info["bonds"] = [[atom_head_id, atom_tail_id, '1', r0, k]]
                #print(atom_head_type, atom_tail_type)
                # 添加angles1
                angles1 = []
                
                #print(bond_info[mol_id_1][1])
                for left_atom in bond_info[mol_id_1][1]:  # 和尾原子链接的原子id和type
                    left_atom_id = str(left_atom[0])
                    left_atom_type = left_atom[1]
                    #angle is left_atom_id-atom_head_id-atom_tail_id
                    # 插入寻找type关系的函数，返回theta0和k
                    angles_p = bond_param['angles']
                    found_param = False
                    for param in angles_p:
                        angle_type1 = left_atom_type + '-' + atom_head_type + '-' + atom_tail_type
                        angle_type2 = atom_head_type + '-' + atom_tail_type + '-' + left_atom_type
                        if param[0] == angle_type1 or param[0] == angle_type2:
                            k = float(param[1]) * 4.184 * 2  # kcal/mol/rad^2 to kJ/mol/rad^2
                            theta0 = float(param[2])
                            found_param = True
                    if found_param == False:
                        print("!!!! No parameter %s or %s was found. Please check the bonged.dat file." % (angle_type1, angle_type2))
                        sys.exit()
                    angles1.append([left_atom_id, atom_head_id, atom_tail_id, '1', theta0, k])
                pair_bond_info["angles1"] = angles1
                # 添加angles2
                angles2 = []
                for right_atom in bond_info[mol_id_2][0]:  # 和首原子链接的原子id和type
                    right_atom_id = str(right_atom[0])
                    right_atom_type = right_atom[1]
                    #angle is atom_head_id-atom_tail_id-right_atom_id
                    # 插入寻找type关系的函数，返回thta0和k
                    angles_p = bond_param['angles']
                    found_param = False
                    for param in angles_p:
                        angle_type1 = right_atom_type + '-' + atom_head_type + '-' + atom_tail_type
                        angle_type2 = atom_head_type + '-' + atom_tail_type + '-' + right_atom_type
                        if param[0] == angle_type1 or param[0] == angle_type2:
                            k = float(param[1]) * 4.184 * 2  # kcal/mol/rad^2 to kJ/mol/rad^2
                            theta0 = float(param[2])
                            found_param = True
                    if found_param == False:
                        print("!!!! No parameter %s or %s was found. Please check the bonged.dat file." % (angle_type1, angle_type2))
                        sys.exit()
                    angles2.append([atom_head_id, atom_tail_id, right_atom_id, '1', theta0, k])
                pair_bond_info["angles2"] = angles2
                # 添加dihedrals
                dihedrals = []
                for left_atom in bond_info[mol_id_1][1]:  # 和尾原子链接的原子id和type
                    left_atom_id = str(left_atom[0])
                    left_atom_type = left_atom[1]
                    for right_atom in bond_info[mol_id_2][0]:  # 和首原子链接的原子id和type
                        right_atom_id = str(right_atom[0])
                        right_atom_type = right_atom[1]
                        # dihedrals is left_atom_id-atom_head_id-atom_tail_id-right_atom_id
                        # 插入寻找type关系的函数，返回func   phase     kd      pn
                        dihedrals_p = bond_param['dihedrals']
                        found_param = False
                        for param in dihedrals_p:
                            dihedrals_type1 = right_atom_type + '-' + atom_head_type + '-' + atom_tail_type + '-' + left_atom_type
                            dihedrals_type2 = left_atom_type + '-' + atom_head_type + '-' + atom_tail_type + '-' + right_atom_type
                            if param[0] == dihedrals_type1 or param[0] == dihedrals_type2:
                                func = int(param[1]) 
                                kd = round(float(param[2]) * 4.184, 2) # kcal/mol to kJ/mol
                                phase = float(param[3]) 
                                pn = float(param[4])
                                found_param = True
                        if found_param == False:
                            for param in dihedrals_p:
                                dihedrals_type3 = 'X' + '-' + atom_head_type + '-' + atom_tail_type + '-' + 'X'
                                dihedrals_type4 = 'X' + '-' + atom_tail_type + '-' + atom_head_type + '-' + 'X'
                                if param[0] == dihedrals_type3 or param[0] == dihedrals_type4:
                                    func = int(param[1]) 
                                    kd = round(float(param[2]) * 4.184, 2) # kcal/mol to kJ/mol
                                    phase = float(param[3]) 
                                    pn = float(param[4])
                                    print("!!! The parameters of %s were used instead of those of %s." % (dihedrals_type3, dihedrals_type1))
                                    found_param = True
                        if found_param == False:
                            print("!!!! No parameter %s or %s was found. Please check the bonged.dat file." % (dihedrals_type1, dihedrals_type2))
                            sys.exit()
                        dihedrals.append([left_atom_id, atom_head_id, atom_tail_id, right_atom_id, func, phase, kd, pn])
                pair_bond_info["dihedrals"] = dihedrals
                #print(pair_bond_info)
                all_bond_info[pair] = pair_bond_info

            # 2. 合并origin_itp和add_itp
            CombineITP(mol_atom_num, repeat_num, "tmp", monomer_itp_name, all_bond_info, res_name, new_atoms, terminal_atoms, origin_atom_id)
            #sys.exit()
            # 3. 保存为res_name的itp, 刷新内存保存, 原始itp进行备份
            #os.system("mv %s.itp %s.itp.bk" % (res_name,res_name))
            os.system("mv tmp.itp %s.itp" % res_name)
            #print(res_name)

            tmp_n += 1
        atom_total_num = natoms + add_atom_num
        #print(add_atom_num)
        #print(atom_total_num)
        #print(totalmoles)
        #sys.exit()
        NewGroName = "md%s.gro" % iter
        WirteGMXGro(totalmoles, NewGroName, atom_total_num, boxsize)
        print("Running round: %s / %s" % (iter, total_iter))
        RunGmx(top_file, iter, em_mdp, md_mdp)
        #sys.exit()
        if iter == total_iter:
            os.system("cp md%s.gro terminated.gro" % iter)
            os.system("rm md*.gro")
            print("*********************************************************************************************")
            print("The termination has been completed, and the terminated structure file is named terminated.gro")
            print("*********************************************************************************************")
        iter += 1

if __name__ == '__main__':
    #args = sys.argv[1:]
    #options = option_parser(args,options)
    em_mdp = "minim.mdp"
    md_mdp = "md.mdp"
    top_file = "top.top"
    gro_file = "init.gro"
    terminal_atoms = ["A1", "D1"]
    new_atoms = ["CA1", "ND1"]
    link_atoms = [1, 1] # 封端链接原子id
    terminal_gros = ["OH.gro","CH3.gro"]
    terminal_itps = ["OH.itp","CH3.itp"]
    bonded_info = "gaff_bonded.dat"
    iter_bond = 3
    print("Start terminating...")
    main(em_mdp, md_mdp, top_file, gro_file, terminal_atoms, new_atoms, link_atoms, terminal_gros, terminal_itps, bonded_info,iter_bond)