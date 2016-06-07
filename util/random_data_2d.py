#!/usr/bin/python

import random

for i in range(30):
	print random.gauss(0.0,0.1), random.gauss(0.0,0.1)

for i in range(30):
	print random.gauss(0.5,0.1), random.gauss(0.5,0.1)

for i in range(30):
	print random.gauss(1.0,0.1), random.gauss(1.0,0.1)

#for i in range(30):
#	print random.gauss(0.0,0.1), random.gauss(5.0,0.1)
