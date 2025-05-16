import sys

structure="structure"
dir = "../Example - Input Pathway/"

xyzr = open(dir + structure + '.xyzr','w') #generate modified file

pqr_file = open(dir + structure + '.pqr','r')
for line in pqr_file:
    if line[:4] == "ATOM":
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]
        r = line[63:70]
        # linemod = x + ' ' + y + ' ' + z + ' ' + r + "\n"
        linemod = x + ' ' + y + ' ' + z + ' ' + r
        #print(linemod)

        xyzr.write(linemod)


#Save output
pqr_file.close()
xyzr.close()

print("Converted to xyzr")
