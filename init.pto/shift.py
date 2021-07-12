#!/usr/bin/env python3
import sys

def within_radius(radius,lattice,pos):
	dx = lattice[0]*0.5-pos[0]
	dy = lattice[1]*0.5-pos[1]
	dz = lattice[2]*0.5-pos[2]
	return dx*dx+dy*dy < radius*radius
	

num_atoms = 0
lattice = "" 

filename = sys.argv[1]
shift = dict()
shift["Ti"] = float(sys.argv[2])
shift["Pb"] = float(sys.argv[3])
shift["O"] = float(sys.argv[4])
radius = 10.0

fout = open('out-'+filename,'w')
with open(filename,'r') as fin:
	num_atoms = int(fin.readline().split()[0])
	line = fin.readline()
	lattice = [float(l) for l in line.split()[0:6]]
	fout.write(f'{num_atoms}\n')
	fout.write(line)
	for line in fin:
		data = line.split()
		elem = data[0]
		pos = [float(x) for x in data[1:4]]

		if pos[0] < 0.5*lattice[0]:
		#if within_radius(radius,lattice,pos):
			if elem == "Ti": pos[2]+=shift["Ti"]
			if elem == "Pb": pos[2]+=shift["Pb"]
			if elem == "O": pos[2]+=shift["O"]

		fout.write(f'{elem} {pos[0]} {pos[1]} {pos[2]}\n')
