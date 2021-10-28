#!/usr/bin/env python3

import shutil, os, argparse, sys, stat, platform
sys.path.append("scripts/pyUtils")
sys.path.append("scripts/setUpScripts")
from utils import Utils
from genFuncs import genHelper
def main():
    name = "njhseq"
    #libs = "bamtools:develop,bibcpp:develop,armadillo:8.200.0"
    libs = "TwoBit:develop,bamtools:develop,boost_math:1_75_0"
    args = genHelper.parseNjhConfigureArgs()
    if Utils.isMac():
        macv, _, _ = platform.mac_ver()
        macv = float('.'.join(macv.split('.')[:2]))
        if macv < 10.12:
            if args.CC and "gcc" in args.CC[0]:
                pass
            else:
                libs = libs + ",sharedMutex:develop"
    cmd = genHelper.mkConfigCmd(name, libs, sys.argv)
    Utils.run(cmd)

main()
