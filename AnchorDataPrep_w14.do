***************************************************************************
***                     PRELIMINARIES                                   ***
***************************************************************************

clear all
set more off		// tells Stata not to pause for --more-- messages
set maxvar 10000	// increases maximal number of variables

global inpath	 "K:\Pairfam\Release14-0\Data\Stata\" // directory of original data


***************************************************************************
***                 EXTRACTING THE VARIABLES                            ***
***************************************************************************

cd $inpath 
use "anchor1.dta"
label language en	

keep id wave sample cohort demodiff intd intm inty east dobm_gen doby_gen age sex_gen nkids cob homosex relstat marstat reldur cohabdur lfs isced2 val1i2 sd30 sd31 pa11 nkids k1type nmar mardur
 										
append using "anchor1_DD.dta", keep(id wave sample cohort demodiff intd intm inty east dobm_gen doby_gen age sex_gen nkids cob homosex relstat marstat reldur cohabdur lfs isced2 val1i2 sd30 sd31 pa11 nkids k1type nmar mardur)

append using "anchor2.dta", keep(id wave sample cohort demodiff intd intm inty east dobm_gen doby_gen age sex_gen nkids cob homosex relstat marstat reldur cohabdur lfs isced2 pa11 nkids k1type nmar mardur)

append using "anchor3.dta", keep(id wave sample cohort demodiff intd intm inty east dobm_gen doby_gen age sex_gen nkids cob homosex relstat marstat reldur cohabdur lfs isced2 val1i2 pa11 nkids k1type nmar mardur)

append using "anchor4.dta", keep(id wave sample cohort demodiff intd intm inty east dobm_gen doby_gen age sex_gen nkids cob homosex relstat marstat reldur cohabdur lfs isced2 pa11 nkids k1type nmar mardur)

* generate combined date of interview
gen intdat=ym(inty, intm)
lab var intdat "Date of interview"
	
gen dateint = ((inty - 1900)*12 + intm - 1)
lab var dateint "Date of interview in months"

save "anchor_w1-4.dta", replace