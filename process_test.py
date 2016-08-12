import argparse, time
from pycompss.api.task import task
from pycompss.api.parameter import *

#class process_test:
#    
#    def __init__(self):
#        self.ready = True
    
@task(x = IN)
def main(x):
    print time.time(), x

y = range(1)

#pt = process_test()
map(main, y)
