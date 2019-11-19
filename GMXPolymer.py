#!/usr/bin/python3
# -*- coding: utf-8 -*
# By Jianchuan Liu  ICCAS  2019-10-15

import linecache
import os
import sys
import math
import numpy as np
import re
import time
import datetime

##############################################################
version="1.0"
authors=["Jianchuan Liu, ICCAS"]
# This is a simple and versatily option class that allows easy
# definition and parsing of options.
class Option:
    def __init__(self, func=str, num=1, default=None, description=""):
        self.func = func
        self.num = num
        self.value = default
        self.description = description

    def __nonzero__(self):
        if self.func == bool:
            return self.value != False
        return bool(self.value)

    def __str__(self):
        return self.value and str(self.value) or ""

    def setvalue(self, v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [self.func(i) for i in v]
###
options = [
    ("-g",        Option(str,                      1,     None, "Input file (GRO)")),
    ("-t",        Option(str,                      1,     None, "Input topology file (TOP)")),
    ("-e",        Option(str,                      1,     None, "Input minim mdp file (MDP)")),
    ("-m",        Option(str,                      1,     None, "Input md mdp file (MDP)")),
    ("-n",        Option(int,                      1,     None, "Input estimate the number of bonds")),
    ("-c",        Option(float,                    1,     None, "Input cut off radius of bond (nm)")),
    ("-i",        Option(str,                      1,     None, "Whether the bond is in the molecule？ Yes (Y) or No (N)")),
    ("-b",        Option(str,                      1,     None, "Input bonds link file (DAT)")),
    ("-h",        Option(bool,                     0,    False, "Display this help")),
    ("-ver",      Option(str,                      1,     None, "Display version")),
    ]
######
def help():
    """Print help text and list of options and end the program."""
    import sys
    print("Usage:\n"
          "python GMXPolymer.py -g init.gro -t topol.top -e minim.mdp "
          "-m md.mdp -n 2 -c 0.6 -i Y -b Bond.dat")
    for item in options:
        if type(item) == str:
            print (item)
    for item in options:
        if type(item) != str:
            print ("%10s  %s"%(item[0],item[1].description))
    print()
    sys.exit()
#
def option_parser(args,options):
    option2 = {}
    option1 = options
    # Check whether there is a request for help
    if len(args) == 0:
        print("Please enter parameters:")
        help()
    ##
    if '-h' in args or '--help' in args:
        help()
    if '-ver' in args or '--version' in args:
        print("                      :-) GMXPloymer, VERSION %s (-:\n"
              "                               2019-10-24\n"
              "                           GMXPloymer is written by:\n"
              "                            %s" % (version,authors[0]))
        sys.exit()
    # Convert the option list to a dictionary, discarding all comments
    options = dict([i for i in options if not type(i) == str])
    #'''
    options['Version'] = version
    options['Arguments'] = args[:]
    while args:
        ar = args.pop(0)
        if len(args) > 0:
            option2[ar] = args[0]
        try:
            options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])
        except Exception as infoerror:
            print("Error in user input:\n"
                  "Invalid command-line\n"
                  "    Unknown command-line: %s" % infoerror)
            sys.exit()

    # This information we would like to print to some files, so let's put it in our information class
    #'''
    for info in option1:
        if info[0] not in options['Arguments'] and info[0] != "-h" and info[0] != "-ver":
            print("!!! Please enter the %s parameters ！！！" % info[0])
            sys.exit()
    #'''
    return option2
    ################
##############################################################
def WroInfo(LineNum):
    '''
    print wrong information
    '''
    print("!"*46)
    print("!!!!! The Gromacs program running error  !!!!!")
    print("!!!!!  The program has been terminated   !!!!!")
    print("!!!!!    All output files are invalid    !!!!!")
    print("!!!!! Please check the code, line %s    !!!!!" % LineNum)
    print("!"*46)
    sys.exit()
    return
#
def WroInfo2(Filename):
    '''
    print wrong information
    '''
    print("!"*55)
    print('This file cannot be found or format error: {:s}'.format(Filename))
    print("!"*55)
    sys.exit()
    return
#
def GmxJudge(CodeLineNum,Number):
    '''
    Judege Whether to execute successfully for Gmx
    '''
    if (os.path.exists('md%s.gro' % Number) != True) or (os.path.getsize('md%s.gro' % Number) == 0):
        WroInfo(CodeLineNum-1)
    '''
    else:
        print("!"*54)
        print("!!!!! The Gromacs program running successfully   !!!!!")
        print("!"*54)
    '''
#
def RunGmx(TopFile,Number,MinimFile,MdFile):
    '''
    Run Gmx program
    '''
    print("Gromacs program is runing...")
    # call gromacs program
    os.system('rm -rf  run.log ./#* *.edr *.cpt *.trr *.tpr md*.gro em*.gro >& /dev/null')
    os.system('gmx grompp -f %s -c tmp%s.gro -p %s -o em%s.tpr -maxwarn 2 >& /dev/null' % (MinimFile,Number-1,TopFile,Number))
    os.system('gmx mdrun -deffnm em%s >& /dev/null' % Number)
    os.system('gmx grompp -f %s -c em%s.gro -p %s -o md%s.tpr -maxwarn 2 >& /dev/null' % (MdFile,Number,TopFile,Number))
    os.system('gmx mdrun -deffnm md%s >& /dev/null' % Number)
    os.system('echo 0 | gmx trjconv -f md%s.gro -s md%s.tpr -o md%s.gro -pbc whole >& /dev/null' % (Number,Number,Number))
    os.system('rm  ./#*  *.edr *.trr *.tpr *.cpt *.xtc em*.gro >& /dev/null')
#
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
def WirteGMXGro(totalmoles, NewGroName, natoms):
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
    for info1 in totalmoles:
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
    NewFile.flush()
    NewFile.close()
    return resid, natom
#
def AddGMXGro(totalmoles, NewGroName, residinput, resname, natomsinput,
              AAtomIndex, DAtomIndex, AAtomType, DAtomType):
    '''
    Add totalmoles: list
    Add  gmx gro file.
    '''
    NewFile = open(NewGroName, 'a')
    n = 1
    m = 1 + natomsinput
    #read total molecle
    resid = n + residinput
    # read single molecle
    for info2 in totalmoles:
        # read single atomic information
        atomtype = info2[1]
        x = info2[2]
        y = info2[3]
        z = info2[4]
        natom = (m - 1) % 99999 + 1
        #print(m,natom,AAtomIndex,DAtomIndex)
        '''
        if natom  == int(AAtomIndex) + m:
            NewFile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %
                          (resid, resname, AAtomType, natom, x, y, z))
        if natom == int(DAtomIndex) + m:
            NewFile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %
                          (resid, resname, DAtomType, natom, x, y, z))
        else:
            NewFile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %
                          (resid, resname, atomtype, natom, x, y, z))
        '''
        NewFile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %
                      (resid, resname, atomtype, natom, x, y, z))
        #'''
        m += 1
    n += 1
    NewFile.flush()
    NewFile.close()
#
def GetDist(AGroLine, DGroLine):
    '''
    Calculate the any two link atom distance
    '''
    xA = float(AGroLine[2])
    yA = float(AGroLine[3])
    zA = float(AGroLine[4])
    xD = float(DGroLine[2])
    yD = float(DGroLine[3])
    zD = float(DGroLine[4])
    ADDist = round(math.sqrt((xA - xD) ** 2 + (yA - yD) ** 2 + (zA - zD) ** 2), 5)
    return ADDist
#
def getThereDimensionListIndex(ThereList, value):
    '''
    get list index of molecule and index of atom
    '''
    MolIndexList = []
    AtomIndexList = []
    for info1 in ThereList:
        for info2 in info1:
            if str(value) == str(info2[1]):
                MolIndex = ThereList.index(info1)
                AtomIndex = info1.index(info2)
                MolIndexList.append(MolIndex)
                AtomIndexList.append(AtomIndex)
    return MolIndexList, AtomIndexList
#
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
        info = re.search(KeyWord, LineData)
        if info == None:
            pass
        else:
            LineList.append(n)
        n += 1
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
    file.close()
    return FileDetail
#
def ReplaceItpFile(filename, keyword, NewType, NewAtom):
    '''
    change the itp file in [ atoms ] information
    :param filename:
    :param keyword:
    :param NewType:
    :param NewAtom:
    :return: change the file
    '''
    LineBegin = GetLineNum(filename, keyword)[0]
    ItpDetail = GetFile(filename)
    if len(LineBegin) != 0:
        AtomLine = ItpDetail[LineBegin[0] - 1]
        AtomLineList = []
        for info in AtomLine.split():
            AtomLineList.append(info)
        AtomLineList[1] = NewType
        AtomLineList[4] = NewAtom
        Empty = ""
        for info in AtomLineList:
            if AtomLineList.index(info) == 0:
                Empty += '%6d' % int(info)
            elif AtomLineList.index(info) == 1:
                Empty += '%5s' % info
            elif AtomLineList.index(info) == 2:
                Empty += '%6d' % int(info)
            elif AtomLineList.index(info) == 3:
                Empty += '%6s' % info
            elif AtomLineList.index(info) == 4:
                Empty += '%6s' % info
            elif AtomLineList.index(info) == 5:
                Empty += '%4d' % int(info)
            elif AtomLineList.index(info) == 6 or AtomLineList.index(info) == 7:
                Empty += '%13.5f' % float(info)
            else:
                Empty += '%10s' % info
        ItpDetail[LineBegin[0] - 1] = Empty + "\n"
        NewItpFile = open(filename, 'w+')
        for info in ItpDetail:
            NewItpFile.write(info)
        NewItpFile.flush()
        NewItpFile.close()
    else:
        pass
#
def InsertInfo(FilName, LineNum, Info="", Position="low"):
    '''
    chang file, add some information in some line
    :param FilName:
    :param LineNum:
    :param Info:
    :param Position:
    :return:
    '''
    if str(Position).lower() == "up":
        LineNum -= 1
        FileDetail = GetFile(FilName)
        FileDetail.insert(LineNum,Info+"\n")
    elif str(Position).lower() == "low":
        FileDetail = GetFile(FilName)
        FileDetail.insert(LineNum,Info+"\n")
    else:
        pass
    filename = open(FilName,'w+')
    for info in FileDetail:
        filename.write(info)
    filename.flush()
    filename.close()
#
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
    return FragmentList,TwoDimensionList
#
def ReadBondLink(BondInforFile):
    linedata = GetItpFrag(BondInforFile, "bonds")[0]
    LinkAD = []
    for line in linedata:
        LinkAD1 = []
        LinkD1 = []
        LinkAatom = str(line.split()[0])
        LinkAtype = str(line.split()[1])
        LinkAnew = str(line.split()[2])
        LinkDatom = str(line.split()[3])
        LinkDtype = str(line.split()[4])
        LinkDnew = str(line.split()[5])
        BondLen = str(line.split()[6])
        BondForce = str(line.split()[7])
        BondFunct = str(line.split()[8])
        LinkAD1.append(LinkAatom)
        LinkAD1.append(LinkAtype)
        LinkAD1.append(LinkAnew)
        LinkAD1.append(LinkDatom)
        LinkAD1.append(LinkDtype)
        LinkAD1.append(LinkDnew)
        LinkAD1.append(BondLen)
        LinkAD1.append(BondForce)
        LinkAD1.append(BondFunct)
        LinkAD.append(LinkAD1)
    return LinkAD
#
def ReadAngleLink(AngleInforFile):
    linedata = GetItpFrag(AngleInforFile, "angles")[0]
    AngleAD = []
    for line in linedata:
        AngleAD1 = []
        AtomType = str(line.split()[0])
        AtomName1 = str(line.split()[1])
        AtomName2 = str(line.split()[2])
        AngleNum = str(line.split()[3])
        ANgleForce = str(line.split()[4])
        Funct = str(line.split()[5])
        AngleAD1.append(AtomType)
        AngleAD1.append(AtomName1)
        AngleAD1.append(AtomName2)
        AngleAD1.append(AngleNum)
        AngleAD1.append(ANgleForce)
        AngleAD1.append(Funct)
        AngleAD.append(AngleAD1)
    return AngleAD
#
def ReadDiheLink(DiheInforFile):
    linedata = GetItpFrag(DiheInforFile, "dihedrals")[0]
    DiheAD = []
    for line in linedata:
        DiheAD1 = []
        AtomType1 = str(line.split()[0])
        AtomName1 = str(line.split()[1])
        AtomName2 = str(line.split()[2])
        AtomType2 = str(line.split()[3])
        Funct = str(line.split()[4])
        C0 = str(line.split()[5])
        C1 = str(line.split()[6])
        C2 = str(line.split()[7])
        C3 = str(line.split()[8])
        C4 = str(line.split()[9])
        C5 = str(line.split()[10])
        DiheAD1.append(AtomType1)
        DiheAD1.append(AtomName1)
        DiheAD1.append(AtomName2)
        DiheAD1.append(AtomType2)
        DiheAD1.append(Funct)
        DiheAD1.append(C0)
        DiheAD1.append(C1)
        DiheAD1.append(C2)
        DiheAD1.append(C3)
        DiheAD1.append(C4)
        DiheAD1.append(C5)
        DiheAD.append(DiheAD1)
    return DiheAD
#
def AllTopfragDict(FileName,TopParentList):
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
#
def WriteNewTop(NewTopName,TopParent,FragDict,AddNum=0):
    '''
    write the new top file accoding link A itp and link D itp
    :param NewTopName:
    :param TopParent:
    :param FragDict:
    :param AddNum:
    :return:
    '''
    FinelName = open(NewTopName,'a')
    if TopParent == "atoms" :
        AtomList = FragDict[TopParent]
        for info in AtomList:
            if (';' in info) == True:
                if info.index(";") == 5:
                    FinelName.write("%6d%5s%6d%6s%6s%6s\n" %
                                    (int(info[0])+AddNum,info[1],int(info[2]),info[3],info[4]),";")
                elif info.index(";") == 6:
                    FinelName.write("%6d%5s%6d%6s%6s%5d%6s\n" %
                                    (int(info[0])+AddNum,info[1],int(info[2]),info[3],info[4],
                                     int(info[5])+AddNum,";"))
                elif info.index(";") == 7:
                    FinelName.write("%6d%5s%6d%6s%6s%5d%13.5f%6s\n" %
                                    (int(info[0])+AddNum,info[1],int(info[2]),info[3],info[4],
                                     int(info[5])+AddNum,float(info[6]),";"))
                elif info.index(";") == 8:
                    FinelName.write("%6d%5s%6d%6s%6s%5d%13.5f%13.5f%6s\n" %
                                    (int(info[0])+AddNum,info[1],int(info[2]),info[3],info[4],
                                     int(info[5])+AddNum,float(info[6]),float(info[7]),";"))
                else:
                    print("The [ %s ] of TmpA.itp or TmpD.itp files format error..., "
                          "please check!!!" % TopParent)
                    sys.exit()
            else:
                if len(info) == 5:
                    FinelName.write("%6d%5s%6d%6s%6s%6s\n" %
                                    (int(info[0])+AddNum,info[1],int(info[2]),info[3],info[4]),";")
                elif len(info) == 6:
                    FinelName.write("%6d%5s%6d%6s%6s%5d%6s\n" %
                                    (int(info[0])+AddNum,info[1],int(info[2]),info[3],info[4],
                                     int(info[5])+AddNum,";"))
                elif len(info) == 7:
                    FinelName.write("%6d%5s%6d%6s%6s%5d%13.5f%6s\n" %
                                    (int(info[0])+AddNum,info[1],int(info[2]),info[3],info[4],
                                     int(info[5])+AddNum,float(info[6]),";"))
                elif len(info) == 8:
                    FinelName.write("%6d%5s%6d%6s%6s%5d%13.5f%13.5f%6s\n" %
                                    (int(info[0])+AddNum,info[1],int(info[2]),info[3],info[4],
                                     int(info[5])+AddNum,float(info[6]),float(info[7]),";"))
                else:
                    print("The [ %s ] of TmpA.itp or TmpD.itp files format error..., "
                          "please check!!!" % TopParent)
                    sys.exit()
    if TopParent == "bonds" :
        AtomList = FragDict[TopParent]
        for info in AtomList:
            if (';' in info) == True:
                # bond type 1 2 4 6
                if info.index(";") == 5:
                    FinelName.write("%6d%7d%4d%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2]),
                                     float(info[3]),float(info[4]),";"))
                # bond type 3
                elif info.index(";") == 6:
                    FinelName.write("%6d%7d%4d%14.4e%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2]),
                                     float(info[3]),float(info[4]),float(info[5]),";"))
                # bond no parameter
                elif info.index(";") == 3:
                    FinelName.write("%6d%7d%4d%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2]),";"))
                else:
                    print("The [ %s ] of TmpA.itp or TmpD.itp files format error..., "
                          "please check!!!" % TopParent)
                    sys.exit()
            else:
                # bond type 1 2 4 6
                if len(info) == 5:
                    FinelName.write("%6d%7d%4d%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2]),
                                     float(info[3]),float(info[4]),";"))
                # bond type 3
                elif len(info) == 6:
                    FinelName.write("%6d%7d%4d%14.4e%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2]),
                                     float(info[3]),float(info[4]),float(info[5]),";"))
                # bond no parameter
                elif len(info) == 3:
                    FinelName.write("%6d%7d%4d%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2]),";"))
                else:
                    print("The [ %s ] of TmpA.itp or TmpD.itp files format error..., "
                          "please check!!!" % TopParent)
                    sys.exit()
    if TopParent == "pairs" :
        AtomList = FragDict[TopParent]
        for info in AtomList:
            if (';' in info) == True:
                if info.index(";") == 3:
                    FinelName.write("%6d%7d%7d%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2]),";"))
                else:
                    print("The [ %s ] of TmpA.itp or TmpD.itp files format error..., "
                          "please check!!!" % TopParent)
                    sys.exit()
            else:
                if len(info) == 3:
                    FinelName.write("%6d%7d%7d%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2]),";"))
                else:
                    print("The [ %s ] of TmpA.itp or TmpD.itp files format error..., "
                          "please check!!!" % TopParent)
                    sys.exit()
    if TopParent == "angles" :
        AtomList = FragDict[TopParent]
        for info in AtomList:
            if (';' in info) == True:
                # angle type 1 2 3
                if info.index(";") == 6:
                    FinelName.write("%6d%7d%7d%7d%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3]),float(info[4]),float(info[5]),";"))
                # angle type 4
                elif info.index(";") == 7:
                    FinelName.write("%6d%7d%7d%7d%14.4e%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3]),float(info[4]),float(info[5]),float(info[6]),";"))
                # angle type 5
                elif info.index(";") == 8:
                    FinelName.write("%6d%7d%7d%7d%14.4e%14.4e%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3]),float(info[4]),float(info[5]),float(info[6]),float(info[7]),";"))
                # angle no parameter
                elif info.index(";") == 4:
                    FinelName.write("%6d%7d%7d%7d%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,
                                     int(info[2])+AddNum,int(info[3]),";"))
                else:
                    print("The [ %s ] of TmpA.itp or TmpD.itp files format error..., "
                          "please check!!!" % TopParent)
                    sys.exit()
            else:
                # angle type 1 2 3
                if len(info) == 6:
                    FinelName.write("%6d%7d%7d%7d%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3]),float(info[4]),float(info[5]),";"))
                # angle type 4
                elif len(info) == 7:
                    FinelName.write("%6d%7d%7d%7d%14.4e%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3]),float(info[4]),float(info[5]),float(info[6]),";"))
                # angle type 5
                elif len(info) == 8:
                    FinelName.write("%6d%7d%7d%7d%14.4e%14.4e%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3]),float(info[4]),float(info[5]),float(info[6]),float(info[7]),";"))
                # angle no parameter
                elif len(info) == 4:
                    FinelName.write("%6d%7d%7d%7d%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,
                                     int(info[2])+AddNum,int(info[3]),";"))
                else:
                    print("The [ %s ] of TmpA.itp or TmpD.itp files format error..., "
                          "please check!!!" % TopParent)
                    sys.exit()
    if TopParent == "dihedrals" :
        AtomList = FragDict[TopParent]
        for info in AtomList:
            if (';' in info) == True:
                # dihedrals type 3
                if info.index(";") == 11:
                    FinelName.write("%6d%7d%7d%7d%7d%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3])+AddNum,int(info[4]),float(info[5]),float(info[6]),
                                     float(info[7]),float(info[8]),float(info[9]),float(info[10]),";"))
                # dihedrals no parameter
                elif info.index(";") == 5:
                    FinelName.write("%6d%7d%7d%7d%7d%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,
                                     int(info[2])+AddNum,int(info[3])+AddNum,int(info[4]),";"))
                # dihedrals type 1
                elif info.index(";") == 8:
                    FinelName.write("%6d%7d%7d%7d%7d%9.2f%10.5f%4d%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3])+AddNum,int(info[4]),float(info[5]),float(info[6]),int(info[7]),";"))
                # dihedrals type 2
                elif info.index(";") == 7:
                    FinelName.write("%6d%7d%7d%7d%7d%9.2f%10.5f%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3])+AddNum,int(info[4]),float(info[5]),float(info[6]),";"))
                else:
                    print("The [ %s ] of TmpA.itp or TmpD.itp files format error..., "
                          "please check!!!" % TopParent)
                    sys.exit()
            else:
                # dihedrals type 3
                if len(info) == 11:
                    FinelName.write("%6d%7d%7d%7d%7d%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3])+AddNum,int(info[4]),float(info[5]),float(info[6]),
                                     float(info[7]),float(info[8]),float(info[9]),float(info[10]),";"))
                # dihedrals no parameter
                elif len(info) == 5:
                    FinelName.write("%6d%7d%7d%7d%7d%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,
                                     int(info[2])+AddNum,int(info[3])+AddNum,int(info[4]),";"))
                # dihedrals type 1
                elif len(info) == 8:
                    FinelName.write("%6d%7d%7d%7d%7d%9.2f%10.5f%4d%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3])+AddNum,int(info[4]),float(info[5]),float(info[6]),int(info[7]),";"))
                # dihedrals type 2
                elif len(info) == 7:
                    FinelName.write("%6d%7d%7d%7d%7d%9.2f%10.5f%6s\n" %
                                    (int(info[0])+AddNum,int(info[1])+AddNum,int(info[2])+AddNum,
                                     int(info[3])+AddNum,int(info[4]),float(info[5]),float(info[6]),";"))
                else:
                    print("The [ %s ] of TmpA.itp or TmpD.itp files format error..., "
                          "please check!!!" % TopParent)
                    sys.exit()
    FinelName.flush()
    FinelName.close()
#
def GetThreeDimListLine(threedimlist,number):
    '''

    :param threedimlist:
    :param number:
    :return:
    '''
    n = 1
    k = 0
    for info in threedimlist:
        if k == number - 1:
            break
        k += 1
        for info1 in info:
            n += 1
    return n
#
def GetMolName(TopFile,AtomNum):
    '''

    :param TopFile:
    :param AtomNum:
    :return:
    '''
    topnum = GetItpFrag(TopFile, "molecules")[1]
    num = 0
    for info in topnum:
        molname = info[0]
        molnum = int(info[1])
        #print(molname,molnum)
        atomtotalnum = int(len((GetItpFrag(molname +".itp", "atoms")[1]))) * molnum
        num += atomtotalnum
        if AtomNum <= num:
            MolName = molname
            break
    return MolName
#
def GetNeighAtom(ItpFile,Aatom,Datom):
    '''

    :param TopFile:
    :param Aatom:
    :param Datom:
    :return:
    '''
    ANeigh = {}
    DNeigh = {}
    AtomList = GetItpFrag(ItpFile, "atoms")[1]
    BondList = GetItpFrag(ItpFile, "bonds")[1]
    for info in BondList:
        Line1 = info[0]
        Line2 = info[1]
        #'''
        # get A atom neighbor
        if int(Aatom) == int(Line1) and int(Datom) != int(Line2):
            ANeigh[Line2] = AtomList[int(Line2)-1][1]
        if int(Aatom) == int(Line2) and int(Datom) != int(Line1):
            ANeigh[Line1] = AtomList[int(Line1)-1][1]
        # get D atom neighbor
        if int(Datom) == int(Line1) and int(Aatom) != int(Line2):
            DNeigh[Line2] = AtomList[int(Line2)-1][1]
        if int(Datom) == int(Line2) and int(Aatom) != int(Line1):
            DNeigh[Line1] = AtomList[int(Line1)-1][1]
        #'''
    return ANeigh, DNeigh



###################################################################################################

def main(options):
    ##############################################################
    # Target total number of bond  in entire system
    TotalBondNum = int(options["-n"]) #2  # integer
    # Bond information file
    BondInforFile = options["-b"] #"Bond.dat"
    # Cut off radius (nm)
    CutOff = float(options["-c"]) #0.6  # nm
    # Whether the bond is in the molecule？ Yes (Y) or No (N)
    BondInMol = options["-i"]#"N"  # Y or N
    # init gmx file
    GroFile = options["-g"]  #"init.gro"
    TopFile = options["-t"]  #"topol.top"
    MinimFile = options["-e"] #"minim.mdp"
    MdFile = options["-m"] #"md.mdp"
    #####################################################
    ##############################################################
    # gmx run max time
    MaxRunGmxTime = 20
    # all top file ParentList
    TopParentList = ["atoms", "bonds", "pairs", "angles", "dihedrals"]
    ##############################################################
    ###################################################################################################
    # initialization
    print("Program begin Runing....")
    # Check if the files exists
    if (os.path.exists(GroFile) != True) or (os.path.getsize(GroFile) == 0):
        WroInfo2(GroFile)
    if (os.path.exists(TopFile) != True) or (os.path.getsize(TopFile) == 0):
        WroInfo2(TopFile)
    if (os.path.exists(MinimFile) != True) or (os.path.getsize(MinimFile) == 0):
        WroInfo2(MinimFile)
    if (os.path.exists(MdFile) != True) or (os.path.getsize(MdFile) == 0):
        WroInfo2(MdFile)
    if (os.path.exists(BondInforFile) != True) or (os.path.getsize(BondInforFile) == 0):
        WroInfo2(BondInforFile)
    ## delete the tmp file
    os.system("rm tmp*.gro >& /dev/null")
    os.system("cp %s tmp0.gro >& /dev/null" % GroFile)
    os.system("cp %s BQ1.top >& /dev/null" % TopFile)
    os.system("rm -rf BondSteep* >& /dev/null")
    # record the begin times
    StartTime = time.time()
    ###
    TotalBondNumBegin = 1
    while TotalBondNumBegin <= TotalBondNum:
        NewGroName = "tmp%s.gro" % TotalBondNumBegin
        filename = "md%s.gro" % TotalBondNumBegin
        TopFile = "BQ%s.top" % TotalBondNumBegin
        NewResname = "BQ%s" % TotalBondNumBegin
        MaxRunGmxTimeBegin = 1
        while MaxRunGmxTimeBegin <= MaxRunGmxTime:
            # Run Gmx first step to pre-equilibrium system
            RunGmx(TopFile,TotalBondNumBegin,MinimFile,MdFile)
            CodeLineNum = sys._getframe().f_lineno
            GmxJudge(CodeLineNum,TotalBondNumBegin)
            # Create a new gro file and new itp file and change the topol.top file
            totalmoles = ReadGMXGro(filename)[2]
            natoms = ReadGMXGro(filename)[0]
            boxsize = ReadGMXGro(filename)[1]
            LinkAD = ReadBondLink(BondInforFile)
            MinDist = float('inf')
            for info1 in LinkAD:
                Aatom = info1[0]
                Datom = info1[3]
                AMolListIndex = getThereDimensionListIndex(totalmoles, Aatom)[0]
                AAtomListIndex = getThereDimensionListIndex(totalmoles, Aatom)[1]
                DMolListIndex = getThereDimensionListIndex(totalmoles, Datom)[0]
                DAtomListIndex = getThereDimensionListIndex(totalmoles, Datom)[1]
                for info2 in AMolListIndex:
                    AtomIndex = AMolListIndex.index(info2)
                    info3 = AAtomListIndex[AtomIndex]
                    AGroLine = totalmoles[info2][info3]
                    for info4 in DMolListIndex:
                        DtomIndex = DMolListIndex.index(info4)
                        info5 = DAtomListIndex[DtomIndex]
                        DGroLine = totalmoles[info4][info5]
                        ADdist = GetDist(AGroLine, DGroLine)
                        if BondInMol.upper() == "Y":
                            if float(ADdist) < MinDist:
                                MinDist = ADdist
                                AMolIndex = info2
                                DMolIndex = info4
                                Aatomname = info1[0]
                                Atype = info1[1]
                                AIndex = info3
                                Anew = info1[2]
                                Datomname = info1[3]
                                Dtype = info1[4]
                                DIndex = info5
                                Dnew = info1[5]
                                BondLen = float(info1[6])
                                BondForce = float(info1[7])
                                Funct = int(info1[8])
                        elif BondInMol.upper() == "N":
                            if float(ADdist) < MinDist and int(info2) != int(info4):
                                MinDist = ADdist
                                AMolIndex = info2
                                DMolIndex = info4
                                Aatomname = info1[0]
                                Atype = info1[1]
                                AIndex = info3
                                Anew = info1[2]
                                Datomname = info1[3]
                                Dtype = info1[4]
                                DIndex = info5
                                Dnew = info1[5]
                                BondLen = float(info1[6])
                                BondForce = float(info1[7])
                                Funct = int(info1[8])
                        else:
                            print("The input is invalid: BondInMol")
                            sys.exit()
            AResname = totalmoles[AMolIndex][0][0]
            DResname = totalmoles[DMolIndex][0][0]
            Anum = GetThreeDimListLine(totalmoles, AMolIndex + 1)
            Dnum = GetThreeDimListLine(totalmoles, DMolIndex + 1)
            Aname = GetMolName(TopFile, Anum)
            Dname = GetMolName(TopFile, Dnum)
            AMol = totalmoles[AMolIndex]
            DMol = totalmoles[DMolIndex]
            ###
            if MinDist <= CutOff:
                print("New bond length is: %s, (%s%s %s %s -> %s%s %s %s), %s bonds" %
                      (MinDist, AMolIndex + 1, AResname, Aatomname, AIndex + 1,
                       DMolIndex + 1, DResname, Datomname, DIndex + 1, TotalBondNumBegin))
                if BondInMol.upper() == "Y" and AMolIndex == DMolIndex:
                    ADMol = AMol
                    ADMolIndexList = []
                    ADMolIndexList.append(AMolIndex)
                    ADMolIndexList.append(DMolIndex)
                    Newtotalmoles = totalmoles
                    Newtotalmoles.pop(max(ADMolIndexList))
                    resid, natom = WirteGMXGro(Newtotalmoles, NewGroName, natoms)
                    AddGMXGro(ADMol, NewGroName, resid, NewResname, natom, AIndex, DIndex, Anew, Dnew)
                    #'''
                    # create new molecule itp file
                    os.system('cp %s.itp %s.itp' % (Aname, NewResname))
                    os.system("sed -i 's/%s/%s/' %s.itp" % (Aname, NewResname, NewResname))
                    # replace itp file of A atom to new type
                    keyword = r'\s'+str(AIndex + 1)+'\s.*\s'+str(NewResname)+'\s.*\s'+str(Aatomname)+'\s'
                    ReplaceItpFile(str(NewResname) + ".itp", keyword, Atype, Anew)
                    # replace itp file of D atom to new type
                    keyword = r'\s'+str(DIndex + 1)+'\s.*\s'+str(NewResname)+'\s.*\s'+str(Datomname)+'\s'
                    ReplaceItpFile(str(NewResname) + ".itp", keyword, Dtype, Dnew)
                    # add the bond information in the itp file
                    LineNum = GetLineNum(str(NewResname) + ".itp", "bonds")[0][0]
                    LineNum += 1
                    Info = "%6d%7d%4d%14.4e%14.4e ; usr define bond" % \
                           (AIndex + 1,DIndex + 1,Funct,BondLen,BondForce)
                    InsertInfo(str(NewResname) + ".itp",LineNum,Info,"low")
                    #'''
                    # add the anlge information in the itp file
                    LineNum = GetLineNum(str(NewResname) + ".itp", "angles")[0][0]
                    AngleAD = ReadAngleLink(BondInforFile)
                    Aneigh = GetNeighAtom(str(NewResname) + ".itp", AIndex + 1, DIndex + 1 + AMolAtomNum)[0]
                    Dneigh = GetNeighAtom(str(NewResname) + ".itp", AIndex + 1, DIndex + 1 + AMolAtomNum)[1]
                    nnn = 1
                    for info6 in Aneigh:
                        atomtype = Aneigh[info6]
                        for info8 in AngleAD:
                            if atomtype == info8[0] and Anew == info8[1] and Dnew == info8[2]:
                                Info = "%6d%7d%7d%7d%14.4e%14.4e ; usr define angle" % \
                                       (int(info6), AIndex + 1,DIndex + 1 + AMolAtomNum,
                                        int(info8[5]),float(info8[3]),float(info8[4]))
                                InsertInfo(str(NewResname) + ".itp",LineNum+1,Info,"low")
                                nnn += 1
                    for info7 in Dneigh:
                        atomtype = Dneigh[info7]
                        for info8 in AngleAD:
                            if atomtype == info8[0] and Dnew == info8[1] and Anew == info8[2]:
                                Info = "%6d%7d%7d%7d%14.4e%14.4e ; usr define angle" % \
                                       (int(info7), DIndex + 1 + AMolAtomNum, AIndex + 1,
                                        int(info8[5]),float(info8[3]),float(info8[4]))
                                InsertInfo(str(NewResname) + ".itp",LineNum+1,Info,"low")
                                nnn += 1
                    # add the dihedral information in the itp file
                    #Aneigh, Dneigh = GetNeighAtom(str(NewResname) + ".itp", AIndex + 1, DIndex + 1 + LineNum)
                    LineNum = GetLineNum(str(NewResname) + ".itp", "dihedrals")[0][0]
                    DiheAD = ReadDiheLink(BondInforFile)
                    for info6 in Aneigh:
                        atomtype1 = Aneigh[info6]
                        for info7 in Dneigh:
                            atomtype2 = Dneigh[info7]
                            for info8 in DiheAD:
                                if atomtype1 == info8[0] and Anew == info8[1] and Dnew == info8[2] and atomtype2 == info8[3]:
                                    Info = "%6d%7d%7d%7d%7d%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e ; usr define dihedral" % \
                                           (int(info6), AIndex + 1,DIndex + 1 + AMolAtomNum,int(info7),
                                            int(info8[4]), float(info8[5]), float(info8[6]),float(info8[7]),
                                            float(info8[8]), float(info8[9]), float(info8[10]))
                                    InsertInfo(str(NewResname) + ".itp", LineNum + nnn, Info, "low")
                                if atomtype1 == info8[3] and Anew == info8[2] and Dnew == info8[1] and atomtype2 == info8[0]:
                                    Info = "%6d%7d%7d%7d%7d%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e ; usr define dihedral" % \
                                           (int(info7), DIndex + 1 + AMolAtomNum, AIndex + 1, int(info6),
                                            int(info8[4]), float(info8[5]), float(info8[6]),float(info8[7]),
                                            float(info8[8]), float(info8[9]), float(info8[10]))
                                    InsertInfo(str(NewResname) + ".itp", LineNum + nnn, Info, "low")
                    #'''
                    # chenge the topol.top file [ molecules ] infromation
                    LineNum,LineTotal = GetLineNum(TopFile, "molecules")
                    LineNum1 = LineNum[0] + 1
                    LineNum2 = LineNum[0] + 1
                    LineNum3 = LineNum[0] + 1
                    AllResList = []
                    while LineNum1 <= LineTotal:
                        ResList = []
                        LineData = str(linecache.getline(TopFile, LineNum1))
                        if re.match(r'\s?;\s.*', LineData) == None and len(LineData.strip()) > 1:
                            Res = LineData.split()[0]
                            Mol = LineData.split()[1]
                            ResList.append(Res)
                            ResList.append(Mol)
                            AllResList.append(ResList)
                        LineNum1 += 1
                    os.system("sed -i '%s,%sd' %s" % (LineNum2,LineTotal,TopFile))
                    for info in AllResList:
                        MolName = info[0]
                        Mole = info[1]
                        if MolName == Aname or MolName == Dname:
                            Mole = int(Mole) - 1
                            os.system("echo ' '%s '      ' %s  >> %s" % (MolName,Mole,TopFile))
                        else:
                            os.system("echo ' '%s '      ' %s  >> %s" % (MolName,Mole,TopFile))
                    os.system("echo ' '%s '      ' %s  >> %s" % (NewResname,1,TopFile))
                    # add the inclue itp in the topol.top
                    LineNum = GetLineNum(TopFile, "system")[0][0]
                    Info = '#include "./%s.itp"' % NewResname
                    InsertInfo(TopFile,LineNum,Info,"up")
                    #'''
                else:
                    ADMol = AMol + DMol
                    ADMolIndexList = []
                    ADMolIndexList.append(AMolIndex)
                    ADMolIndexList.append(DMolIndex)
                    Newtotalmoles = totalmoles
                    Newtotalmoles.pop(max(ADMolIndexList))
                    Newtotalmoles.pop(min(ADMolIndexList))
                    resid, natom = WirteGMXGro(Newtotalmoles, NewGroName, natoms)
                    AddGMXGro(ADMol, NewGroName, resid, NewResname, natom, AIndex, DIndex, Anew, Dnew)
                    #'''
                    # create new molecule itp file
                    os.system('cp %s.itp tmpA%s.itp' % (Aname,TotalBondNumBegin))
                    os.system('cp %s.itp tmpD%s.itp' % (Dname,TotalBondNumBegin))
                    # change the mol name
                    os.system("sed -i 's/%s/%s/' %s%s.itp" % (Aname, NewResname, "tmpA" ,TotalBondNumBegin))
                    os.system("sed -i 's/%s/%s/' %s%s.itp" % (Dname, NewResname, "tmpD" ,TotalBondNumBegin))
                    # Merge the two itps
                    AMolAtomNum = len(GetItpFrag(str("tmpA%s.itp" % TotalBondNumBegin), "atoms")[0])
                    AFragDict = AllTopfragDict("tmpA"+str(TotalBondNumBegin)+".itp",TopParentList)
                    DFragDict = AllTopfragDict("tmpD"+str(TotalBondNumBegin)+".itp",TopParentList)
                    NewItpFile = open(NewResname + ".itp" , 'w+')
                    NewItpFile.write("; User Define itp file\n")
                    NewItpFile.write("\n")
                    NewItpFile.write("\n")
                    NewItpFile.write("[ moleculetype ]\n")
                    NewItpFile.write(";name            nrexcl\n")
                    NewItpFile.write("%4s%15d\n" % (NewResname, 3))
                    NewItpFile.write("\n")
                    NewItpFile.flush()
                    NewItpFile.close()
                    for info in TopParentList:
                        os.system("echo [ %s ]  >> %s" % (info,NewResname + ".itp"))
                        # write A mol information
                        WriteNewTop(NewResname + ".itp",info,AFragDict,0)
                        # write D mol information
                        WriteNewTop(NewResname + ".itp",info,DFragDict,AMolAtomNum)
                    # replace itp file of A atom to new type
                    keyword = r'\s'+str(AIndex + 1)+'\s.*\s'+str(NewResname)+'\s.*\s'+str(Aatomname)+'\s'
                    ReplaceItpFile(str(NewResname) + ".itp", keyword, Atype, Anew)
                    # replace itp file of D atom to new type
                    keyword = r'\s'+str(DIndex + 1 + AMolAtomNum)+'\s.*\s'+str(NewResname)+'\s.*\s'+str(Datomname)+'\s'
                    ReplaceItpFile(str(NewResname) + ".itp", keyword, Dtype, Dnew)
                    # add the bond information in the itp file
                    LineNum = GetLineNum(str(NewResname) + ".itp", "bonds")[0][0]
                    Info = "%6d%7d%4d%14.4e%14.4e ; usr define bond" % \
                           (AIndex + 1,DIndex + 1 + AMolAtomNum,Funct,BondLen,BondForce)
                    InsertInfo(str(NewResname) + ".itp",LineNum,Info,"low")
                    #'''
                    # add the anlge information in the itp file
                    LineNum = GetLineNum(str(NewResname) + ".itp", "angles")[0][0]
                    AngleAD = ReadAngleLink(BondInforFile)
                    Aneigh = GetNeighAtom(str(NewResname) + ".itp", AIndex + 1, DIndex + 1 + AMolAtomNum)[0]
                    Dneigh = GetNeighAtom(str(NewResname) + ".itp", AIndex + 1, DIndex + 1 + AMolAtomNum)[1]
                    nnn = 1
                    for info6 in Aneigh:
                        atomtype = Aneigh[info6]
                        for info8 in AngleAD:
                            if atomtype == info8[0] and Anew == info8[1] and Dnew == info8[2]:
                                Info = "%6d%7d%7d%7d%14.4e%14.4e ; usr define angle" % \
                                       (int(info6), AIndex + 1,DIndex + 1 + AMolAtomNum,
                                        int(info8[5]),float(info8[3]),float(info8[4]))
                                InsertInfo(str(NewResname) + ".itp",LineNum+1,Info,"low")
                                nnn += 1
                    for info7 in Dneigh:
                        atomtype = Dneigh[info7]
                        for info8 in AngleAD:
                            if atomtype == info8[0] and Dnew == info8[1] and Anew == info8[2]:
                                Info = "%6d%7d%7d%7d%14.4e%14.4e ; usr define angle" % \
                                       (int(info7), DIndex + 1 + AMolAtomNum, AIndex + 1,
                                        int(info8[5]),float(info8[3]),float(info8[4]))
                                InsertInfo(str(NewResname) + ".itp",LineNum+1,Info,"low")
                                nnn += 1
                    # add the dihedral information in the itp file
                    #Aneigh, Dneigh = GetNeighAtom(str(NewResname) + ".itp", AIndex + 1, DIndex + 1 + LineNum)
                    LineNum = GetLineNum(str(NewResname) + ".itp", "dihedrals")[0][0]
                    DiheAD = ReadDiheLink(BondInforFile)
                    for info6 in Aneigh:
                        atomtype1 = Aneigh[info6]
                        for info7 in Dneigh:
                            atomtype2 = Dneigh[info7]
                            for info8 in DiheAD:
                                if atomtype1 == info8[0] and Anew == info8[1] and Dnew == info8[2] and atomtype2 == info8[3]:
                                    Info = "%6d%7d%7d%7d%7d%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e ; usr define dihedral" % \
                                           (int(info6), AIndex + 1,DIndex + 1 + AMolAtomNum,int(info7),
                                            int(info8[4]), float(info8[5]), float(info8[6]),float(info8[7]),
                                            float(info8[8]), float(info8[9]), float(info8[10]))
                                    InsertInfo(str(NewResname) + ".itp", LineNum + nnn, Info, "low")
                                if atomtype1 == info8[3] and Anew == info8[2] and Dnew == info8[1] and atomtype2 == info8[0]:
                                    Info = "%6d%7d%7d%7d%7d%14.4e%14.4e%14.4e%14.4e%14.4e%14.4e ; usr define dihedral" % \
                                           (int(info7), DIndex + 1 + AMolAtomNum, AIndex + 1, int(info6),
                                            int(info8[4]), float(info8[5]), float(info8[6]),float(info8[7]),
                                            float(info8[8]), float(info8[9]), float(info8[10]))
                                    InsertInfo(str(NewResname) + ".itp", LineNum + nnn, Info, "low")
                    #'''
                    # chenge the topol.top file [ molecules ] infromation
                    LineNum,LineTotal = GetLineNum(TopFile, "molecules")
                    LineNum1 = LineNum[0] + 1
                    LineNum2 = LineNum[0] + 1
                    LineNum3 = LineNum[0] + 1
                    AllResList = []
                    while LineNum1 <= LineTotal:
                        ResList = []
                        LineData = str(linecache.getline(TopFile, LineNum1))
                        if re.match(r'\s?;\s.*', LineData) == None and len(LineData.strip()) > 1:
                            Res = LineData.split()[0]
                            Mol = LineData.split()[1]
                            ResList.append(Res)
                            ResList.append(Mol)
                            AllResList.append(ResList)
                        LineNum1 += 1
                    os.system("sed -i '%s,%sd' %s" % (LineNum2,LineTotal,TopFile))
                    for info in AllResList:
                        MolName = info[0]
                        Mole = info[1]
                        if MolName == Aname or MolName == Dname:
                            Mole = int(Mole) - 1
                            os.system("echo ' '%s '      ' %s  >> %s" % (MolName,Mole,TopFile))
                        else:
                            os.system("echo ' '%s '      ' %s  >> %s" % (MolName,Mole,TopFile))
                    os.system("echo ' '%s '      ' %s  >> %s" % (NewResname,1,TopFile))
                    # add the inclue itp in the topol.top
                    LineNum = GetLineNum(TopFile, "system")[0][0]
                    Info = '#include "./%s.itp"' % NewResname
                    InsertInfo(TopFile,LineNum,Info,"up")
                #'''
                # wirte the box size in the new gro file
                NewGroFile = open(NewGroName, 'a')
                NewGroFile.write("%10.5f%10.5f%10.5f\n" % (boxsize[0], boxsize[1], boxsize[2]))
                NewGroFile.close()
                # back off infromation
                os.system("mkdir BondSteep-%s >& /dev/null" % TotalBondNumBegin)
                os.system("cp *.gro *.top *.itp ./BondSteep-%s" % TotalBondNumBegin)
                # delete the tmp file
                os.system("rm tmp*.itp tmp*.itp >& /dev/null")
                break
            #
            else:
                os.system("rm tmp%s.gro tmp%s.gro >& /dev/null" % (TotalBondNumBegin - 1, TotalBondNumBegin))
                os.system("cp md%s.gro tmp%s.gro >& /dev/null " % (TotalBondNumBegin,TotalBondNumBegin-1))
                print("No new bond were product. continue runing gromacs, %s/%s Time" % (MaxRunGmxTimeBegin, MaxRunGmxTime))
            MaxRunGmxTimeBegin += 1
        # Determine the number of times the Gromacs program has run
        if MaxRunGmxTimeBegin > MaxRunGmxTime:
            print("The Gromacs program has reached its maximum number of runs")
            print("There are formed %s bond, %s bonds are still unformed" %
                  (TotalBondNumBegin-1,TotalBondNum-TotalBondNumBegin+1))
            os.system("mkdir BondSteepFailed >& /dev/null")
            os.system("cp *.gro *.top *.itp ./BondSteepFailed >& /dev/null")
            os.system("cp ./BondSteepFailed/BQ%s.top ./BondSteepFailed/Failed.top >& /dev/null" % (TotalBondNumBegin))
            os.system("cp ./BondSteepFailed/md%s.gro ./BondSteepFailed/Failed.gro >& /dev/null" % (TotalBondNumBegin))
            os.system('rm -rf *.log ./#* *.xtc *.edr *.cpt *.trr *.tpr mdout.mdp '
                      'md*.gro em*.gro tmp*.gro BQ* >& /dev/null')
            print("!!!The all Failed file is in the 'BondSteepFailed' folder!!!")
            EndTime = time.time()
            TotalTime = EndTime-StartTime
            HMS = datetime.timedelta(seconds=TotalTime)
            print('!!! it costs %s !!!' % HMS)
            sys.exit()
        os.system("mv BQ%s.top BQ%s.top" % (TotalBondNumBegin,TotalBondNumBegin + 1))
        TotalBondNumBegin += 1
    ###################################################################################################
    ##
    if TotalBondNumBegin > TotalBondNum:
        os.system("mkdir BondSteepSuccess >& /dev/null")
        os.system('gmx grompp -f %s -c tmp%s.gro -p BQ%s.top -o Success.tpr -maxwarn 2 >& /dev/null' %
                  (MinimFile,TotalBondNumBegin-1,TotalBondNumBegin))
        os.system('gmx mdrun -deffnm Success >& /dev/null')
        os.system('gmx grompp -f %s -c Success.gro -p BQ%s.top -o Success.tpr -maxwarn 2 >& /dev/null' %
                  (MdFile,TotalBondNumBegin))
        os.system('gmx mdrun -deffnm Success >& /dev/null')
        os.system("cp *.gro *.top *.itp ./BondSteepSuccess")
        os.system("cp ./BondSteepSuccess/BQ%s.top ./BondSteepSuccess/Success.top" % (TotalBondNumBegin))
        print("!!!                Program completion              !!!")
        print("There are formed %s bond, %s bonds are still unformed" %
              (TotalBondNumBegin-1,TotalBondNum-TotalBondNumBegin+1))
        print("!!!The all successful file is in the 'BondSteepSuccess' folder!!!")
        EndTime = time.time()
        TotalTime = EndTime-StartTime
        HMS = datetime.timedelta(seconds=TotalTime)
        print('!!! it costs %s !!!' % HMS)
        os.system('rm -rf BondSteep-* *.log ./#* *.xtc *.edr *.cpt *.trr *.tpr mdout.mdp '
                  'md*.gro em*.gro tmp*.gro Success*.gro BQ* >& /dev/null')
        sys.exit()

##########################################################################################################
if __name__ == '__main__':
    args = sys.argv[1:]
    options = option_parser(args,options)
    main(options)
# The End

