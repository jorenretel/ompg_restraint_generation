

import os
import subprocess
import inspect

def get_git_hash():
    caller_file = inspect.stack()[1][1]
    dirname = os.path.dirname(os.path.realpath(caller_file))
    git_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=dirname)

    if len(git_hash) == 8:
        return git_hash[:-1]
    else:
        return ''