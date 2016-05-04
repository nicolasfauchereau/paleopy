def save_progress(path=None, step=None, value=0):
    import os
    import json
    from collections import OrderedDict as od
    """
    save a progress indicator (`step`) in
    a JSON file
    """
    progress = od()
    progress['step'] = step
    progress['percentage'] = value
    with open(os.path.join(path, 'output.json'), 'w') as f:
        json.dump(progress, f)
