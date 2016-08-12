#!/usr/bin/env python

import shutil, os, argparse, sys, stat
sys.path.append("scripts/pyUtils")
sys.path.append("scripts/setUpScripts")
from utils import Utils
from genFuncs import genHelper
def main():
    name = "bibseq"
    libs = "bamtools:v2.4.0,bibcpp:release/v2.4.0,armadillo:6.200.3,TwoBit:release/v2.0.4,hts:1.3.1"
    args = genHelper.parseNjhConfigureArgs()
    if Utils.isMac():
        if args.CC and "gcc" in args.CC[0]:
            pass
        else:
            libs = libs + ",sharedMutex:release/v0.6"
    cmd = genHelper.mkConfigCmd(name, libs, sys.argv, "-lcurl")
    Utils.run(cmd)
main()
