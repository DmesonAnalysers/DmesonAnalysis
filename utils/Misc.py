import math
import sys

def exec_time(start, stop):
    '''
    Prints the exection time in hh:mm:ss. 'time' library is required.

    Usage
    ----------
    start = time.time()

    ...
    some code
    ...

    exec_time(start, time.time())
    '''

    hours = math.floor((stop - start)/3600)
    minutes = math.floor((stop - start - hours*3600)/60)
    seconds = math.floor((stop - start - hours*3600 - minutes*60))

    print(f'Execution time: {hours}:{minutes}:{seconds}')

def Debug(msg = '', printMessageInfo=True):
    '''
    Prints the file name and the line in the code in which this function is called with a message
    '''
    if printMessageInfo:
        print(f"debug file {sys.argv[0]} line {sys._getframe().f_back.f_lineno}: {msg}")