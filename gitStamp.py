'''This module is used to find the current git hash in the
   directory (of the python file from which this function
   is called). It is handy to 'stamp' a created dataset
   with this hash, so it is clear in the future with which
   version of the script, the dataset was created. In order
   for this to work, the latest changes should be committed.

'''

import os
import subprocess
import inspect

def get_git_hash():
    '''Get git hash from the directory of the file is located
       the caller file is located in.

    '''

    caller_file = inspect.stack()[1][1]
    dirname = os.path.dirname(os.path.realpath(caller_file))
    git_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'],
                                       cwd=dirname)

    if len(git_hash) == 8:
        return git_hash[:-1]
    else:
        return ''
