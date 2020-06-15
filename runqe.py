#!/usr/bin/python3
import os

def runpw(inp, out):
	return os.system("mpirun -n 4 pw.x <" + inp + "> " + out)
