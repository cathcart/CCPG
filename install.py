#!/usr/bin/env python
import shutil
#this is an install script for the cross code pseudopotential generator code

def edit_make_depend():
  location="make.depend"
  current=open(location).read().strip().split("\n")
  dict={}#this will hold a list of the current files and their dependencies
  for i in current:
    t=[x.strip() for x in i.split(":")]
    if t[0] not in dict.keys():
      dict[t[0]]=[t[1]]
    else:
      dict[t[0]]=dict[t[0]]+[t[1]]

  #now check to see if we've already updated the file
  install_files=["siesta_pass.a","gsl_interface_f.o","gsl_interface_c.o"]
  shutil.copy2(location,location+".backup")
  file_out=open(location,'a')
  for file in install_files:
    if file not in dict["gener_pseudo.o"]:
      file_out.write("gener_pseudo.o : %s\n"%file)
      print "...Added %s dependency to make.depend." %file
    else:
      print "...%s dependency already in make.depend." %file
  print "...make.depend updated."
  file_out.close()

def edit_Makefile():
  location="Makefile"
  flags=[False,False]#flags to test if the Makefile has been correctly edited

  current=open(location).read().strip().split("\n")
  dict={}#this will hold a list of the current lines
  for i in current:
    t=[x.strip() for x in i.split()]
    if len(t)>0 and t[0] not in dict.keys():
      dict[t[0]]=t[1::]

  if "LD1LIBS=" in dict.keys():#ok
    flags[0]=True
    print "...LD1LIBS defintion included"
  if "$(LD1LIBS)" in dict["ld1.o"]:#ok
    flags[1]=True
    print "...ld1.o line ok"

  shutil.copy2(location,location+".backup")

  if flags[0] and flags[1]:
    print "...The Makefile is correct"
  else:
    file_out=open(location,'w')
    for line in open(location+".backup").read().strip().split("\n"):
      t=line.split()

      if flags[0] == False:
        if len(t)>0 and t[0]=="TLDEPS=":
          line+="\n\nLD1LIBS= -lgsl -lgslcblas siesta_pass.a"
          print "...Added LD1LIBS defintion"

      if flags[1] == False:
        if len(t)>0 and t[0]=="ld1.o":
          line +=" $(LD1LIBS)"
          print "...Added LD1LIBS to complile line"

      file_out.write(line+"\n")

def edit_gener_pseudo():
  location="gener_pseudo.f90"
  
  current=open(location).read()

  #first include
  insert_text="\n  use siesta_pass, only: siesta_output"
  if current.find(insert_text) >= 0:
    print "...use siesta_pass module already included"
  else:
    a=current.find("use ")
    b=current.find("\n",a)
    current=current[:b]+insert_text+current[b:]
    print "...use siesta_pass module added"

  #second include
  insert_text="elseif (pseudotype == 4) then\n     write(stdout, &\n          '(/,5x,21(''-''),'' Generating SIESTA pseudopotential '',21(''-''),/)')\n     pseudotype =1\n  "
  if current.find(insert_text) >= 0:
    print "...new pseudotype definition already included"
  else:
    a=current.find("else\n")
    current=current[:a]+insert_text+current[a:]
    print "...new pseudotype definition added"

  #third include
  insert_text="\n  if (pseudotype .eq. 1) then\n   call siesta_output(grid)\n  endif"
  if current.find(insert_text) >= 0:
    print "...siesta_output call already included"
  else:
    a=current.rfind("endif\n")
    b=current.find("\n",a)
    current=current[:b]+insert_text+current[b:]
    print "...siesta_output call added"

  print "...gener_pseudo.f90 updated"

  #write output
  shutil.copy2(location,location+".backup")
  file_out=open(location,'w')
  file_out.write(current)

def edit_ld1_readin():
  location="ld1_readin.f90"

  current=open(location).read()
  insert_text="if (pseudotype < 1.or.pseudotype > 4)"  
  if current.find(insert_text) >= 0:
    print "...ld1_readin.f90 already corrected"
  else:
    a=current.find("f (pseudotype < 1.or.pseudotype > 3)")
    current=current[:a]+insert_text+current[a+len(insert_text):]
    print "...ld1_readin.f90 corrected"

  #write output
  shutil.copy2(location,location+".backup")
  file_out=open(location,'w')
  file_out.write(current)

if __name__=="__main__":
  print "make.depend"
  #edit_make_depend()
  print "Makefile"
  #edit_Makefile()
  print "gener_pseudo.f90"
  #edit_gener_pseudo()
  print "Installation complete"
 
  edit_ld1_readin()
