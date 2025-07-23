---------------------------------------------------------------
| "Analysis of time-to-event data with two time scales.       |
|  An application to transitions out of cohabitation"         | 
| Angela Carollo, Hein Putter, Paul Eilers and Jutta Gampe    |
---------------------------------------------------------------

July 22, 2025

This is a repository for the code to reproduce the analyses
presented in the paper "Analysis of time-to-event data with two 
time scales. An application to transitions out of cohabitation"

This repository consists of many files, mostly written in R and
one file written in Stata.

The data used in this paper come from the "Panel Analysis of 
Intimate Relationships and Family Dynamics" (pairfam), 
waves 1-11. 
The data cannot be shared publicly, but can be applied for at
https://www.pairfam.de/en/data/data-access/

The code is extensively commented, but if any question arises
please write at: carollo@demogr.mpg.de

Note: Figures 1-3 and Figures A1-A2 are illustrative, and 
      present no substantive data findings. We do not provide 
      scripts to reproduce them.
---------------------------------------------------------------
* Data used in this project:
 - Release 11: 'biopart.dta'
 - Release 11: 'anchor1.dta', 'anchor1_DD.dta', 'anchor2.dta',
	       'anchor3.dta', 'anchor4.dta', 'anchor5.dta',
               'anchor4.dta', 'anchor5.dta', 'anchor6.dta',
               'anchor6.dta', 'anchor7.dta', 'anchor8.dta',
               'anchor9.dta', 'anchor10.dta', 'anchor11.dta'	

---------------------------------------------------------------
* R packages used in this project:
 - readstata13 (to read .dta files into R)
 - Hmisc       (for the function 'describe()')
 - data.table  (for easy manipulation of the datasets and vars)
 - TwoTimeScales (for the analysis)
 - fields      (for 'image.plot()')
 - viridis     (for color palettes)
 - JOPS        (for bbase, used in VCM calculations)

---------------------------------------------------------------
* source files used in this project:
 - "iwsl_1d.R"
 - "VCM_2TS_fixed.R"

---------------------------------------------------------------
* How to run the code:
1. Download the data and place them in a directory
2. Run the file "Prep_stata.do" to prepare the wave's specific
   data (containing the covariates of interests)
3. Run "1_data_preparation.R" to prepare the main dataset for 
   the analysis. This file cleans the 'biopart' data and also
   incorporates the covariates from waves 1-11.
   The main file for the analysis is now "one_cohabitation18.Rda"

---- MODELS FIT ----
4. The scripts "2.1_tts_marriage.R" and "2.2_tts_separation.R"
   fit the two time scales P-splines models to the events 
   marriage and separation separately, for the four subgroups 
   of Men and Women, East and West Germany - without covariates.
   They produce the models in Tables 3 and 4,
   (third column "full bivariate").
5. The scripts "3.1_add_marriage.R" and "3.2_add_separation.R"
   fit the log-additive (P-GAM) P-splines models to the events
   marriage and separation separately for the four subgroups.
   They produce the models in Tables 3 and 4, first column
   ("log-additive").
6. The scripts "4.1_vcm_marriage.R" and "4.2_vcm_separation.R"
   fit the varying-coefficients models with P-splines to the
   events marriage and separation separately for the four groups.
   They produce the models in Tables 3 and 4, second column
   ("varying-coefficients).

---- FIGURES (main paper) ----
7. The script "5_Figure4.R" reproduces the two panels in Figure 4
8. The script "6_Figure5.R" reproduces Figure 5
9. The script "7_Figure6.R" reproduces Figure 6


