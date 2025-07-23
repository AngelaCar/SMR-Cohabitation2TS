* CHANGE ALL PATHS TO YOUR DATA FOLDER!!!
use "K:\Pairfam\Release11-0\Data\Stata\anchor1.dta" 

keep id wave intd intm inty doby_gen dobm_gen east val1i2 isced isced2 sd30 sd31 										
append using "K:\Pairfam\Release11-0\Data\Stata\anchor1_DD.dta", keep(id wave intd intm inty doby_gen dobm_gen east val1i2 isced isced2 sd30 sd31)
append using "K:\Pairfam\Release11-0\Data\Stata\anchor2.dta", keep(id wave intd intm inty doby_gen dobm_gen east isced isced2)
append using "K:\Pairfam\Release11-0\Data\Stata\anchor3.dta", keep(id wave intd intm inty doby_gen dobm_gen east val1i2 isced isced2)
append using "K:\Pairfam\Release11-0\Data\Stata\anchor4.dta", keep(id wave intd intm inty doby_gen dobm_gen east isced isced2)
append using "K:\Pairfam\Release11-0\Data\Stata\anchor5.dta", keep(id wave intd intm inty doby_gen dobm_gen east val1i2 isced isced2 sd30 sd31)
append using "K:\Pairfam\Release11-0\Data\Stata\anchor6.dta", keep(id wave intd intm inty doby_gen dobm_gen east isced isced2)
append using "K:\Pairfam\Release11-0\Data\Stata\anchor7.dta", keep(id wave intd intm inty doby_gen dobm_gen east val1i2 isced isced2)
append using "K:\Pairfam\Release11-0\Data\Stata\anchor8.dta", keep(id wave intd intm inty doby_gen dobm_gen east isced isced2)
append using "K:\Pairfam\Release11-0\Data\Stata\anchor9.dta", keep(id wave intd intm inty doby_gen dobm_gen east val1i2 isced isced2 sd30 sd31)
append using "K:\Pairfam\Release11-0\Data\Stata\anchor10.dta", keep(id wave intd intm inty doby_gen dobm_gen east isced isced2)
append using "K:\Pairfam\Release11-0\Data\Stata\anchor11.dta", keep(id wave intd intm inty doby_gen dobm_gen east val1i2 isced isced2 sd30 sd31)

* generate combined date of interview
	gen intdat=ym(inty, intm)
	lab var intdat "Date of interview"
	
gen dateint = ((inty - 1900)*12 + intm - 1)
lab var dateint "Date of interview in months"

gen bdate = ((doby_gen - 1900)*12 + dobm_gen -1)
lab var bdate "Date of birth in months"

gen ageint = (dateint - bdate)
lab var ageint "Age at interview in months"

gen ageint_y = ageint/12
lab var ageint_y "Age at interview in years"
 
* change the following path to your data folder!!!!
save "C:\Cloud\PhD project\Applications\P2_CohabitationMarriageSeparation\Analysis_Wave1-11\SMR_analysis_new\Analysis\anchor_all_waves.dta", replace

