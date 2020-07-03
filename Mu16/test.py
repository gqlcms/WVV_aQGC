code_string = ""
for j in (open("EDBR2PKUTree.h", 'r').readlines()):
  for i in (open("hhh.txt", 'r').readlines()):
    var = i.split("->")[1].replace(" ","").replace("\n","")
    if len(var)>0:
       if (var in j) & ("SetBranchAddress" in j):
          j = "\n"
  code_string += j

with open("new.h","w") as f:
     f.write(code_string)
