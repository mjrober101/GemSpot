#initialize varibles x for reading in file and pdb for actual atoms
import os
x = []
pdb=[]
ligand = []
endres=[]
startres=[1]
output = open("corrected.pdb","w")
output2 = open("lig.pdb","w")
forjaws = open("pepzinput1","w")
forjaws2 = open("pepzinput2","w")
allcommands = open ("allcommands.sh","w")

#read in maestro pdb file
with open('/path/to/project/choppedpdb.pdb') as f:
	for line in f:
		if line[:4] == 'ATOM' or line[:6] == "HETATM":
			splitline = [line[:6], line[7:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54],line[54:60],line[60:66],line[66:80]]
			splitline=[elem.strip() for elem in splitline]
			x.append(splitline)
#pull out lines corresponding to actual atoms
def pdbfilter(file,string,out):
	for i in range(0,len(file)):
		y = file[i][0]
		if y == string:
			out.append(file[i])

pdbfilter(x,"ATOM",pdb)
pdbfilter(x,"HETATM",pdb)

an=1
ligatomnum=1
renumber=1
#write out pdb file, inserting TER records in between, renumbering from 1, deleting segid, renaming ligand atom types

for i in range(0,len(pdb)):	
	#segment pdb information into each column
	cat = (pdb[i][0])
	atomnum = int((pdb[i][1]))
	atomname = (pdb[i][2])
	resname = (pdb[i][3])
	chainid = (pdb[i][4])
	resid = int(pdb[i][5])
	xcoor = float(pdb[i][6])
	ycoor = float(pdb[i][7])
	zcoor = float(pdb[i][8])
	occu = float(pdb[i][9])
	bfact = float(pdb[i][10])
	scat = (pdb[i][11])
	#these next 4 are used for chain breaks
	if i < len(pdb)-1:
		a=i+1	
	residf = int(pdb[a][5])
	chainidf = (pdb[a][4])
	#write pdb
	if cat =="ATOM":
		output.write("%-5s %5i  %-4s%3s %1s %3i %11.3f %7.3f %7.3f %5.2f%5.2f %11s\n" % (cat, an, atomname, resname, " ",renumber, xcoor, ycoor, zcoor, occu, bfact, scat))
		an = (an + 1)
		#check for breaks, insert TER when found
		if residf - resid > 1:
			output.write("TER\n")
			#keep track of and write length of chains e.g. 1-20, 21-40. Needed for PEPZ input file
			endres.append(renumber)
			startres.append(renumber + 1)
		#check for changes in segment, inesrt TER when found
		elif chainid != chainidf:
			output.write("TER\n")
			endres.append(renumber)
			startres.append(renumber + 1)
		#iterate up the resnumber when the residue switches starting from 1
		if residf - resid != 0:
			renumber = renumber + 1
	elif cat=="HETATM":
		output.write("%6s%5i  %s%02d %3s %1s %3i %11.3f %7.3f %7.3f %5.2f%5.2f %11s\n" % (cat, an, atomname[0],ligatomnum, resname, " ",renumber, xcoor, ycoor, zcoor, occu, bfact, scat))
		output2.write("%6s%5i  %s%02d %3s %1s %3i %11.3f %7.3f %7.3f %5.2f%5.2f %11s\n" % (cat, an, atomname[0],ligatomnum, resname, " ",1, xcoor, ycoor, zcoor, occu, bfact, scat))
		ligatomnum = ligatomnum + 1
		an = an + 1
endres.append(renumber)

bosspath=os.environ["BOSSdir"]
ligcommand = bosspath+"scripts/xPDBZ lig"
allcommands.write("#!/bin/bash\n")
allcommands.write("%s\n" % (ligcommand))
ligcommand2 = bosspath+"scripts/xSPM lig"
allcommands.write("%s\n" % (ligcommand2))
allcommands.write('grep -B 500 Geometry lig.z >ligx.z\n')
allcommands.write("grep -A 500 'Variable Bonds follow' sum >> ligx.z\n")

#now we write out the lines for the pepz input. This is a fortran program,
#so max line length of 80 must be enforced. So, we will write in blocks of
#four segements. First input generates fully flexible zmat for optimization
#def pepzwriter(destination,reslist1,reslist2):

def makefortranloop(destination,string,reslist1,reslist2):
	j=0
	destination.write("\n%s" % (string))
	for i in range(0,len(reslist1)):
		if j < 4:
			destination.write("%s%i%s%i" % (" ",reslist1[i],"-",reslist2[i]))
			j=j+1
		else:
			j=1
			destination.write("\n%s" % (string))
			destination.write("%s%i%s%i" % (" ",reslist1[i],"-",reslist2[i]))

def makepepzinput(destination,string1,string2,string3,reslist1,reslist2,fixed):
	destination.write("$ title Jaws from glideEM\n")
	destination.write("$ read database oplsaam.db\n")
	destination.write("$ read dihedrals dihedrals.aam\n")
	destination.write("$ read parameter oplsaam.par")
	makefortranloop(destination,string1,reslist1,reslist2)
	makefortranloop(destination,string2,reslist1,reslist2)
	destination.write("\n$ read boss ligx.z\n")
	destination.write("$ read pdb corrected.pdb\n")
	destination.write("$ center")
	makefortranloop(destination,string3,reslist1,reslist2)
        if fixed == 1:
		makefortranloop(destination,"$ set fixed backbone",reslist1,reslist2)
	destination.write("\n$ write zmatrix foropt.z\n")
	destination.write("$ write pdb pepz1.pdb\n")

makepepzinput(forjaws,"$ set override domain","$ set parameter type ALL","$ set variable all",startres,endres,0)
makepepzinput(forjaws2,"$ set override domain","$ set parameter type ALL","$ set variable bonds,angles,improper,unsaturated",startres,endres,1)
#this writes a pepz input that will be used to generate a fixed backbone
#zmat. 
mcpropath=os.environ['MCPROdir']
pepzcommand = mcpropath+"pepz/pepz -i pepzinput1 -o pepzoutput1"
allcommands.write("%s\n" % (pepzcommand))
optcommand = mcpropath +"scripts/xCG92 foropt"
allcommands.write("%s\n" % (optcommand))
allcommands.write('grep -B 9999 Geometry optsum >Jaws.z\n')
allcommands.write('rm foropt.z\n')
allcommands.write('rm pepz1.pdb\n')
pepzcommand2 = mcpropath+"pepz/pepz -i pepzinput2 -o pepzoutput2"
allcommands.write("%s\n" % (pepzcommand2))
allcommands.write("grep -A 99999 'Variable Bonds follow' foropt.z >>Jaws.z\n")
