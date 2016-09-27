"""
Copyright 2016 EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from pycompss.api.parameter import *
from pycompss.api.task import task

def main(x):
    from pycompss.api.api import compss_wait_on
    print "Main process:"
    results = []
    for i in x:
        results.append(print_time(i))
    results = compss_wait_on(results)
    print results

@task(x = IN, returns = int)
def print_time(x):
    import time
    x = time.time()
    return x

if __name__ == "__main__":
    y = range(10)
    main(y)
