#!/usr/bin/env python3
import os
import sys
from ruamel.yaml import YAML

with open(os.path.expanduser(sys.argv[1]), 'r') as f:
    data = YAML().load(f)

data = dict(data)
for sample in data.keys():
    values = dict(data[sample])
    print("sample={}, values={}".format(sample, values))
