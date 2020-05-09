## small script for producing doxygen files from doctest files
import os

for fileName in os.listdir("."):
  if fileName.endswith(".doctest"):
    txtFile = open(fileName,'r')
    doxyFile= open(fileName[:-7]+"doxy",'w')
    for line in txtFile:
      if "#doctest: +SKIP" in line:  line=line.replace('#doctest: +SKIP','')
      if ", doctest=True)" in line:  line=line.replace(', doctest=True)',')')
      if "doctestImage("   in line:
        line=line.replace('>>> doctestImage("','')
        line=line.replace('")\n','')
        line="\\endverbatim\n"+\
              "\\image html "+line+".png width=40%\n"+\
              "\\verbatim"
      elif "doctest"       in line:  continue
      doxyFile.write(line)
    txtFile.close()
    doxyFile.close()
