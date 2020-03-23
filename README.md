# smle - Simulation-based estimation of peer effects (Krauth 2006)

This repository includes the SMLE package of programs for estimating the model in my 2006 Journal of Econometrics
paper "Simulation-based estimation of peer effects."
It includes Windows binaries, Fortran source code, and a Stata 
package.

## Installing the package in Windows

### Stata

The current general release of the Stata package can be obtained by executing the 
Stata command:
```stata
net install smle, from("http://www.sfu.ca/~bkrauth/code")
```
To see how the command works, you can call  `help smle`

The current developmental version can be obtained from this site by executing the Stata command:

```stata
net install smle, from("https://raw.githubusercontent.com/bvkrauth/smle/master/stata/")
```

### Binaries only

While the Stata package is the most user-friendly way to estimate the model, 
it is only a wrapper for the binary executables 
smle.exe and s2.exe. If you do not have access to Stata, you can 
download the binaries from  https://github.com/bvkrauth/smle/releases
and run them directly.

See [doc/smledoc.pdf](doc/smledoc.pdf) for additional information.

## Installing the package in other operating systems

If you are not using Windows, you can download and compile the source code
for the Fortran programs.

The program has been compiled successfully with a variety of systems and 
compilers.

See [doc/smledoc.pdf](doc/smledoc.pdf) for additional information.

## Folder structure

The folder structure of this repository is:

  - src: source code
  - doc: documentation
  - stata: stata package
  - testing: a few test files
 