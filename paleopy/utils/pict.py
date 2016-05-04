import os
from collections import OrderedDict as od
import json

def save_progress(path=None, step=None, value=0):
    progress = od()
    progress['step'] = step
    progress['percentage'] = value
    with open(os.path.join(path, 'output.json'), 'w') as f:
        json.dump(progress, f)

