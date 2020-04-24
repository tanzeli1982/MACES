#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 18:14:58 2020

Modify the specified setting in the namelist file 

@author: Zeli Tan
"""

import sys
import xml.etree.ElementTree as ET
from optparse import OptionParser

# use the OptionParser to parse the command line options
parser = OptionParser()
parser.add_option("-f", "--file", type="string", dest="filename")
parser.add_option("-p", "--parameter", type="string", dest="parameter")
parser.add_option("-v", "--value", type="string", dest="value")

(options, args) = parser.parse_args()

# extract arguments and check
filename = options.filename
parameter = options.parameter
value = options.value

try:
    isfind = False
    # find the node and modify its value
    tree = ET.parse(filename)
    root = tree.getroot()
    findstr = "./group/entry"
    for entry in root.findall(findstr):
        if entry.get('id')==parameter:
            isfind = True
            entry.set('value', value)
    assert isfind, "parameter " + parameter + " is not found in the xml file"
    # write to xml file
    tree.write(filename)
except AssertionError as errstr:
    # print error message and exit the program
    print("xmlchange fails due to that", errstr)
    sys.exit()