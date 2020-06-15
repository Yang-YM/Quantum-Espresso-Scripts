#!/usr/bin/python3
import re

class pwinp:
	inp = []
	parameter = {}
	data = {}
	orderP = []
	orderD = []
	nat = 0
	ntyp = 0
	atom_type = []
	# read pw.x input file and preprocess
	def read(self, name):
		self.inp = []
		self.parameter = {}
		self.data = {}
		self.orderP = []
		self.orderD = []
		self.atom_type = []
		file = open(name)
		temp = file.read()
		file.close()
		temp = temp.replace("'", "").replace("\"", "")
		content = re.split("\n|,", temp)
		i = 0
		start = 0
		end = 0
		while i < len(content):
			if len(content[i]) == 0:
				i += 1
				continue
			if content[i][0] == "&":
				start = i
				content[i].replace(" ", "")
				i += 1
				while content[i][0] != "/":
					i += 1
				i += 1
				end = i
				self.inp.append(content[start:end])
				self.orderP.append(content[start][1:])
				self.parameter[self.orderP[-1]] = self.inp[-1]
			else:
				i += 1
		self.nat = int(self.findP("SYSTEM", "nat")[0][1])
		self.ntyp = int(self.findP("SYSTEM", "ntyp")[0][1])
		i = end
		while i < len(content):
			if len(content[i]) < 8:
				i += 1
				continue
			content[i].replace(" ", "")
			if content[i][0:8] == "ATOMIC_S":
				self.orderD.append("ATOMIC_SPECIES")
				self.inp.append(content[i:i + self.ntyp + 1])
				i += self.ntyp + 1
				self.data[self.orderD[-1]] = self.inp[-1]
			elif content[i][0:8] == "ATOMIC_P":
				self.orderD.append("ATOMIC_POSITIONS")
				self.inp.append(content[i:i + self.nat + 1])
				i += self.nat + 1
				self.data[self.orderD[-1]] = self.inp[-1]
			elif content[i][0:8] == "CELL_PAR":
				self.orderD.append("CELL_PARAMETERS")
				self.inp.append(content[i:i + 4])
				i += 4
				self.data[self.orderD[-1]] = self.inp[-1]
			elif content[i][0:8] == "K_POINTS":
				self.orderD.append("K_POINTS")
				temp = re.split("{|}|\(|\)", content[i])
				if temp[1] == "gamma":
					self.inp.append(content[i:i + 1])
					i += 1
				elif temp[1] == "automatic":
					self.inp.append(content[i:i + 2])
					i += 2
				else:
					self.temp = int(content[i + 1]) + 1
					self.inp.append(content[i:i + temp + 1])
					i += temp + 1
				self.data[self.orderD[-1]] = self.inp[-1]
			elif content[i][0:8] == "ATOMIC_F":
				self.orderD.append("ATOMIC_FORCES")
				self.inp.append(content[i:i + self.nat + 1])
				i += self.nat + 1
				self.data[self.orderD[-1]] = self.inp[-1]
			elif content[i][0:8] == "CONSTRAI":
				self.orderD.append("CONSTRAINTS")
				temp = re.split("{|}|\(|\)", content[i + 1])
				temp = int(temp[0]) + 1
				self.inp.append(content[i:i + temp + 1])
				i += temp + 1
				self.data[self.orderD[-1]] = self.inp[-1]
			else:
				i += 1
		for i in range(1, self.ntyp + 1):
			self.atom_type.append(self.data["ATOMIC_SPECIES"][i].split()[0])
	# write pw.x input file
	def write(self, name):
		contentP = []
		contentD = []
		for i in self.orderP:
			contentP.append("\n".join(self.parameter[i]))
		for i in self.orderD:
			contentD.append("\n".join(self.data[i]))
		file = open(name, 'w')
		file.write("\n".join(contentP))
		file.write("\n\n")
		file.write("\n\n".join(contentD))
		file.close()
	# find all calculation parameters with gived 'tag' in a 'part'
	def findP(self, part, tag):
		p = self.parameter[part]
		r = []
		for i in range(1, len(p) - 1):
			if 0 < p[i].find(tag) < p[i].find("="):
				temp = p[i].split()
				temp[1] = temp[2]
				temp[2] = i
				r.append(temp[0:3])
		return r
	# remove all calculation parameters with gived 'tag' in a 'part'
	def removeP(self, part, tag):
		p = self.parameter[part]
		n = len(p) - 1
		for i in range(1, n):
			if p[n - i].find(tag) != -1:
				del p[n - i]
	# insert calculation parameters in a position of a 'part'
	def insertP(self, part, para, pos = -1):
		p = self.parameter[part]
		flag = 0
		if pos == -1:
			pos = len(p) - 1
		for i in range(0, len(para)):
			flag = 0
			temp = self.findP(part, para[i][0])
			for j in range(0, len(temp)):
				if temp[j][0] == para[i][0]:
					flag = 1
					if temp[j][1] != para[i][1]:
						temp_ = p[temp[j][2]].split("=")
						p[temp[j][2]] = temp_[0] + "= " + para[i][1]
			if flag == 0:
				p.insert(pos, "    " + para[i][0] + " = " + para[i][1])
				pos += 1
