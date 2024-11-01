# Programed in Python3 but run in python2
# -*- coding: utf-8 -*-
"fork!!! from https://github.com/CarlosCHMS/CSU.git"

"""
Created on Wed Mar 18 15:39:49 2020

@author: carlos
"""
#python3 CSU.py ./caseFolder/

import sys
import subprocess as sp
import os

if __name__=="__main__":
    exec_path = './bin/FVM2D'
    if not os.path.exists(exec_path):
        raise RuntimeError('please compile source file!')
    if len(sys.argv) < 1:
        raise ValueError('please choose case!')  
    case_path = sys.argv[1]
    if not os.path.exists(case_path):
        raise ValueError( 'Case path is incorrect or does not exist!')  
    if case_path[-1] != '/':                                        # 用于拼接路径
        case_path += '/'
    print("\nFVM2D - A CFD code for unstructured meshs:\n")    
    os.system("rm  %ssolution.csv" % case_path)                     # 删除之前运行过的结果
    os.system("%s %s" % (exec_path, case_path))                     # CFD计算
    os.system("python3 %sanalisys.py %s" % (case_path, case_path))  # 后处理

