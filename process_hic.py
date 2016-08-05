#!/usr/bin/python
# -*- coding: utf-8 -*-
'''process Hi-C paired end FastQ files'''
import time
from pycompss.api.task import task
from pycompss.api.parameter import *

from fastq2adjacency import fastq2adjacency


