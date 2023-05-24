# import numpy as np
import sys
import os.path

filename = sys.argv[1]
init = int(sys.argv[2])
fin = int(sys.argv[3])
#filename_output = sys.argv[4]

hist, hist_reac, scored, ne, ng = 0,0,0,0,0
mine, ming, maxe, maxg = [] ,[] ,[] ,[]

#output = open("./work/results/"+filename+".phsp", "w")
notfound=0
merged=0
noempty=0
for i in range(init,fin+1,1):
    try:
        if (os.path.exists("./work/run"+str(i)+"/log.err")):
            f = open("./work/run"+str(i)+"/"+filename+".header", "r")
#            f2 = open("./work/run"+str(i)+"/"+filename+".phsp", "r")
            mensaje = f.readlines()
            if len(mensaje) > 0:
#                data = f2.read()
#                output.write(data)
                hist = hist + int(mensaje[2][30:])
                hist_reac = hist_reac + int(mensaje[3][54:])
                scored = scored + int(mensaje[4][27:])
                ne = ne + int(mensaje[18][13:])
                ng = ng + int(mensaje[19][17:])
                mine.append(float(mensaje[21][30:-5]))
                ming.append(float(mensaje[22][32:-5]))

                maxe.append(float(mensaje[24][29:-5]))
                maxg.append(float(mensaje[25][32:-5]))
                noempty+=1
            f.close()
            merged+=1
        else: notfound+=1
    except:
        print(f'Problems processing file run{i}/{filename}')

f = open("./work/results/"+filename+".header", "w")
f.write("TOPAS ASCII Phase Space\n")
f.write("\n")
f.write("Number of Original Histories: "+str(hist)+"\n")
f.write("Number of Original Histories that Reached Phase Space: "+str(hist_reac)+"\n")
f.write("Number of Scored Particles: "+str(scored)+"\n")
f.write("\n")
f.write("Columns of data are as follows:")
f.write("1: Position X [cm]\n")
f.write("2: Position Y [cm]\n")
f.write("3: Position Z [cm]\n")
f.write("4: Direction Cosine X\n")
f.write("5: Direction Cosine Y\n")
f.write("6: Energy [MeV]\n")
f.write("7: Weight\n")
f.write("8: Particle Type (in PDG Format)\n")
f.write("9: Flag to tell if Third Direction Cosine is Negative (1 means true)\n")
f.write("10: Flag to tell if this is the First Scored Particle from this History (1 means true)\n")
f.write("\n")
f.write("Number of e-: "+str(ne)+"\n")
f.write("Number of gamma: "+str(ng)+"\n")
f.write("\n")
f.write("Minimum Kinetic Energy of e-: "+str(min(mine))+" MeV\n")
f.write("Minimum Kinetic Energy of gamma: "+str(min(ming))+" MeV\n")
f.write("\n")
f.write("Maximum Kinetic Energy of e-: "+str(max(maxe))+" MeV\n")
f.write("Maximum Kinetic Energy of gamma: "+str(max(maxg))+" MeV\n")
f.close()

print(str(notfound)+" jobs are still running \n")
print(str(merged)+" jobs have been merged \n")
print(f'{merged-noempty} empty files')
print(f'{noempty} files merged successfuly.')


