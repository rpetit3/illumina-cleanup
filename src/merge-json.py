#! /usr/bin/env python
import json
import os
import sys
from collections import OrderedDict

json_data = OrderedDict()
for arg in sys.argv:
    if arg.endswith(".json"):
        name = os.path.basename(arg).split(".")[0]
        with open(arg) as fh:
            json_data[name] = json.load(fh, object_pairs_hook=OrderedDict)

print(json.dumps(json_data, indent=4))
