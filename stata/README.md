# smle/stata - Stata package for smle`

This folder contains the developmental version of the SMLE Stata package 

Folder contents are:
 - smle.ado: the ado-file for the smle command
 - smle.hlp: the help file for the smle command
 - smle_example.do: a Stata do-file that demonstrates and tests the package
 
You can install the package by hand. Just put the files smle.ado, smle.hlp, and the current 
version of the binary executables smle.exe and s2.exe in either your working directory or 
in your ADOFILEs\s directory.

Alternatively you can install it within Stata by executing the command:
```stata
net install smle, from("https://raw.githubusercontent.com/bvkrauth/smle/master/stata/")
```

 
