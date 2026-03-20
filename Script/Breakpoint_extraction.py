#!/usr/bin/python3
# -*- coding: utf-8 -*
# By Jianchuan Liu  XHU  2026-03-20
import os

def get_integer_from_user():
    """
    持续提示用户输入一个整数，直到输入有效。
    返回用户输入的整数。
    """
    while True:
        user_input = input("输入你想从哪一个断点提取可模拟的文件(整数): ").strip()
        if not user_input:
            print("输入不能为空，请重新输入。")
            continue
        try:
            value = int(user_input)
            if value <= 0:
                print("输入不能是0或负数, 请输入一个正整数。")
                continue
            return value
        except ValueError:
            print("输入无效，请输入一个整数。")
Breakid = get_integer_from_user()
################################################################################################################
print("!!!The all file is in the 'BreakMD_%s' folder!!!" % (int(Breakid)))
print("!!!The box file will named 'init.gro'; the top file will named 'top.top'!!!")
os.system("mkdir BreakMD_%s >& /dev/null" % (int(Breakid)))
os.system("cp -rf BondSteep-%s/md*.gro ./BreakMD_%s/init.gro " % (int(Breakid), int(Breakid)))
os.system("cp -rf BondSteep-%s/B%s.top ./BreakMD_%s/top.top " % (int(Breakid)-1, int(Breakid)-1, int(Breakid)))
os.system("cp -rf BondSteep-%s/*.itp ./BreakMD_%s " % (int(Breakid)-1, int(Breakid)))



