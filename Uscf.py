#!/usr/bin/python3

import sys
import os
import re
import libpwinp as pw
import runqe

# SCF calculation with given alpha
def scf(alpha):
	global inp, Utype, file_name, ext
	para = []
	for i in range(1, len(Utype)):
		if Utype[i] == 1:
			para.append(["Hubbard_alpha(" + str(i) + ")", str(alpha)])
	inp.insertP("SYSTEM", para)
	inp.write(file_name + ext)
	temp = runqe.runpw(file_name + ext, file_name + ".out")
	if temp != 0:
		print("ERROR, please check input and output files")
		sys.exit()

# process output and return on site occupations
def occupations():
	global file_name
	file = open(file_name + ".out")
	out = file.read()
	file.close()
	temp = out.split(" --- enter write_ns ---\n")
	temp0 = temp[2].split("\n --- exit write_ns ---")[0].split("\n")
	temp1 = temp[3].split("\n --- exit write_ns ---")[0].split("\n")
	start = 0
	length = 0
	for i in range(0, len(temp0)):
		if temp0[i][0:5] == "atom ":
			if start == 0:
				start = i
			else:
				length = i - start
				break
	num = (len(temp0) - start) // length
	atom = ["    "] * num
	n0 = [0.0] * num
	n1 = [0.0] * num
	for i in range(0, num):
		atom[i] = int(temp0[start].split()[1])
		n0[i] = float(temp0[start].split()[-1])
		n1[i] = float(temp1[start].split()[-1])
		start += length
	return atom, n0, n1

# program
file_name = sys.argv[1]
ext = file_name[file_name.rfind("."):]
file_name = file_name[:-len(ext)]
inp = pw.pwinp()
inp.read(file_name + ext)
# check input file
inp.insertP("CONTROL", [["calculation", "'scf'"]])
temp = inp.findP("SYSTEM", "Hubbard_U")
if len(temp) == 0:
	print("Please add Hubbard_U parameter")
	sys.exit()
inp.insertP("SYSTEM", [["lda_plus_u", ".true."]])
# obtain Hubbard U calculation information
temp = inp.findP("CONTROL", "prefix")
save = ""
if len(temp) == 0:
	inp.insertP("CONTROL", [["prefix", "'" + file_name + "-U'"]])
	save = file_name + "-U.save"
else:
	save = temp[0][1] + ".save"
type = {}
for i in range(0, inp.ntyp):
	type[inp.atom_type[i]] = i + 1
Utype = [0] * (inp.ntyp + 1)
U = [0.0] * (inp.ntyp + 1)
equal = []
temp = inp.findP("SYSTEM", "Hubbard_U")
for i in temp:
	num = int(re.split("\(|\)", i[0])[1])
	Utype[num] = 1
	U[num] = float(i[1])
display = ""
for i in range(1, len(Utype)):
	if Utype[i] == 1:
		display += str(i) + " (" + inp.atom_type[i - 1] + ") "
print("Calculate the self-consistent Hubbard U")
print("Atomic types: " + display)
print("\nConstrain Hubbard U. If not, please continue.")
print("e.g. \"1 2\": Fix U(1) U(2), \"1=2=3\": Make U(1)=U(2)=U(3)")
temp = input("Atomic types: ").split()
for i in temp:
	if i.find("=") != -1:
		equal.append([])
		temp_ = i.split("=")
		for j in temp_:
			if 0 < int(j) < len(Utype):
				equal[-1].append(int(j))
	elif 0 < int(i) < len(Utype):
		Utype[int(i)] = -1
# Hubbard U iteration
iter = 0
maxdiff = 1
while maxdiff > 0.1:
	iter += 1
	print("\nIteration " + str(iter))
	print("Hubbard U for calculation")
	for i in range(1, len(Utype)):
		if Utype[i] != 0:
			print("%4d"%i + "%5s"%inp.atom_type[i - 1] + " : " + "%2.2f"%U[i] + " eV")
	para = []
	for i in range(1, len(Utype)):
		if Utype[i] == 1:
			para.append(["Hubbard_U(" + str(i) + ")", "%.2f"%U[i]])
	inp.insertP("SYSTEM", para)
	# SCF calculation for alpha = 0
	print("SCF calculation alpha =  0")
	scf(0.00)
	os.system("cp -r " + save + " " + save + "-U")
	# SCF calculation for alpha = 0.05
	inp.insertP("ELECTRONS", [["startingpot", "file"], ["startingwfc", "file"]])
	print("SCF calculation alpha =  0.05")
	scf(0.05)
	atom, a1n0, a1n1 = occupations()
	os.system("rm -rf " + save)
	os.system("mv " + save + "-U " + save)
	# SCF calculation for alpha = -0.05
	print("SCF calculation alpha = -0.05")
	scf(-0.05)
	atom, a2n0, a2n1 = occupations()
	# calculate Hubbard U with linear response
	Ucal = []
	for i in range(0,len(Utype)):
		Ucal.append(0.0)
		Ucal.append([])
	d = inp.data["ATOMIC_POSITIONS"]
	print("linear response Hubbard U")
	print("  atom  type  Hubbard U")
	for i in range(0, len(atom)):
		type_ = d[atom[i]].split()[0]
		U_ = 0.1 / (a1n0[i] - a2n0[i]) - 0.1 / (a1n1[i] - a2n1[i])
		print("%5d"%atom[i] + "%6s"%type_ + "%9s"%"%2.5f"%U_)
		Ucal[2 * type[type_] + 1].append(U_)
	for i in range(1, len(Utype)):
		if len(Ucal[2 * i + 1]) != 0:
			Ucal[2 * i] = sum(Ucal[2 * i + 1]) / len(Ucal[2 * i + 1])
	# results
	print("Calculated Hubbard U")
	for i in range(1, len(Utype)):
		if len(Ucal[2 * i + 1]) != 0:
			print("%4d"%i + "%5s"%inp.atom_type[i - 1] + " : " + "%2.2f"%Ucal[2 * i] + " eV" + \
		 " (max: " + "%2.2f"%max(Ucal[2 * i + 1]) + ", min: "+ "%2.2f"%min(Ucal[2 * i + 1]) + " )")
	# update Hubbard U
	for i in equal:
		Usum = 0
		for j in i:
			Usum += Ucal[2 * j]
		Usum /= len(i)
		for j in i:
			Ucal[2 * j] = Usum
	maxdiff = 0
	for i in range(1, len(Utype)):
		if Utype[i] == 1:
			if maxdiff < abs(U[i] - Ucal[2 * i]):
				maxdiff = abs(U[i] - Ucal[2 * i])
			U[i] += 0.66 * (Ucal[2 * i] - U[i])
# end of iteration
print("\nFinish Hubbard U iteration")
print("Final Hubbard U")
for i in range(1, len(Utype)):
	if Utype[i] != 0:
		print("%4d"%i + "%5s"%inp.atom_type[i - 1] + " : " + "%2.2f"%U[i] + " eV")
print("SCF calculation for final Hubbard U")
para = []
for i in range(1, len(Utype)):
	if Utype[i] == 1:
		para.append(["Hubbard_U(" + str(i) + ")", "%.2f"%U[i]])
inp.insertP("SYSTEM", para)
inp.removeP("SYSTEM", "Hubbard_alpha")
inp.write(file_name + ext)
runqe.runpw(file_name + ext, file_name + ".out")