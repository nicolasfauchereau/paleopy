def save_progress(path=None, step=None, value=0):
    """
    save a progress indicator (`step`) in
    a JSON file
    """
    progress = od()
    progress['step'] = step
    progress['percentage'] = value
    with open(os.path.join(path, 'output.json'), 'w') as f:
        json.dump(progress, f)
