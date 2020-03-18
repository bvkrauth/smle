discard
clear
sysuse census, clear

gen own = (popurban > 3328253)
gen npeers = region
gen peeravg = round(runiform()*region)/region

smle own peeravg npeers pop , replace execute

* bootstrap _b , reps(5): smle own peeravg npeers pop , replace execute

