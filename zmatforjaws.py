#initialize varibles x for reading in zmat
x = []
ligatoms =[]
ligsites = []
output = open("jawsfinal.z","w")
forjaws0 = open ("phase0.par","w")
forjaws1 = open ("phase1.par","w")
forjaws2 = open ("phase2.par","w")
with open('/path/to/project/Jaws.z') as f:
	for line in f:
		x.append(line.split())

#Write new z matrix without TERZ and storing ligand atoms
for i in range(1,len(x)):
	if x[i][0] != "Geometry":
		atomnum = int(x[i][0])
		atomname = (x[i][1])
		atomtype1 = int(x[i][2])
		atomtype2 = int(x[i][3])
		bondedatom = int(x[i][4])
		distance = float(x[i][5])
		angleatom = int(x[i][6])
		angle = float(x[i][7])
		diheatom = int(x[i][8])
		dihe = float(x[i][9])
		resid = (x[i][10])
		resnum = int(x[i][11])
#		if resid != "LIG":
		output.write("%4i%4s %4i %4i %4i %11.6f%4i %11.6f%4i %11.6f %3s %4i\n" % \
		(atomnum, atomname, atomtype1, atomtype2, bondedatom, distance, angleatom, angle, diheatom, dihe, resid, resnum))
		if resid == "LIG":
			ligatoms.append(int(x[i][0]))
	elif x[i][0] == "Geometry": 
		break

#Add in TERZ and CAP atom for centering 
output.write("TERZ\n")
atomnum = atomnum+1
resnum = resnum +1
output.write("%4i%4s %4i %4i %4i %11.6f%4i %11.6f%4i %11.6f %3s %4i\n" % \
		(atomnum, "CAP", -1, -1, bondedatom, 1.0, angleatom, 120.0, diheatom, 180.0, "CAP", resnum))

liglength=float(len(ligatoms)-1)
h=liglength /10
print ligatoms
for i in range(0,10):
	m = (int(round(1+h*i)))
	ligsites.append(ligatoms[m])
