#!/usr/bin/env python

import shutil, os, argparse, sys, stat, platform
sys.path.append("scripts/pyUtils")
sys.path.append("scripts/setUpScripts")
from utils import Utils
from genFuncs import genHelper
def main():
    name = "bibseq"
    libs = "bamtools:release/v2.5.0,bibcpp:release/v2.5.0,armadillo:7.600.1"
    args = genHelper.parseNjhConfigureArgs()
    if Utils.isMac():
        macv, _, _ = platform.mac_ver()
        macv = float('.'.join(macv.split('.')[:2]))
        if macv < 10.12:
            if args.CC and "gcc" in args.CC[0]:
                pass
            else:
                libs = libs + ",sharedMutex:release/v0.7"
    cmd = genHelper.mkConfigCmd(name, libs, sys.argv)
    
    Utils.run(cmd)
    
main()

