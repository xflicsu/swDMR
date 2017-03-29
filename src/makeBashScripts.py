#!/usr/bin/env python
# encoding: utf-8
"""
makeBashScripts.py
"""

import sys
import os


def main():
    tool_map =  {'intersect': 'intersectBed',}
    # create a BASH script for each old tool, mapping to the new CLI command.
    for tool in tool_map:
        new = tool
        old = tool_map[tool]
        
        script = open('bin/'  + old, 'w')
        script.write("#!/bin/sh\n")
        script.write("${0%/*}/bedtools " + new + " \"$@\"\n")
        script.close()

if __name__ == "__main__":
    main()
