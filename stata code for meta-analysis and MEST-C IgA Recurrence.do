
///////////////////////////////////////////////////////////////////////////////
**#  START META_ANALYSIS O CRUDE AND ADJ EST. ALACHKAR/BAGNASCO vs PARK
///////////////////////////////////////////////////////////////////////////////
clear
cd "C:\Documenti\Delsante\Bagnasco IgA"

////////////////////////////////////////////////////////////////////////////////
**# Input Crude/Adj Est. & 95%CI from Table 2 (pp2861) Park AJT 2019;19:2855
////////////////////////////////////////////////////////////////////////////////

preserve
clear
input str8 comp  n d perc hr1 lb1 ub1 pval1 hr2 lb2 ub2 pval2 hr3 lb3 ub3 pval3
"M" 35 77 45.5 2.34 1.52 3.58 0.0001 1.68 1.04 2.71 0.03 1.59 0.98 2.58 0.06
"E" 33 99 33.3 1.66 1.08 2.57 0.02 1.66 1.04 2.63 0.03 1.55 0.96 2.51 0.07
"S" 52 148 35.1 2.38 1.55 3.64 0.0001 1.92 1.16 3.20 0.01 1.81 1.06 3.11 0.03
"T1" 22 62 35.5 2.33 1.39 3.90 0.001 1.67 0.91 3.08 0.10 1.63 0.88 3.03 0.12
"T2" 23 42 54.8 3.70 2.23 6.14 0.0001 2.90 1.65 5.08 0.0001 2.85 1.61 5.01 0.0001
"C" 16 44 36.4 1.76 1.02 3.04 0.04 1.23 0.67 2.23 0.50 1.13 0.61 2.10 0.71
"C2" 4   6 66.7 5.51 1.99 15.25 0.001 5.75 1.84 17.96 0.002 5.06 1.57 16.31 0.007
"ATCMR" 35 146 24.0 1.04 0.68 1.60 0.86 0.98 0.61 1.58 0.94 1.08 0.68 1.72 0.75
"ABMR" 23 78 29.5 1.65 1.02 2.67 0.04 1.53 0.88 2.68 0.13 1.63 0.93 2.87 0.09
end
* list
save park_table2, replace

////////////////////////////////////////////////////////////////////////////////
**# Reshape the Park's Est % 95% CI from Model 1,2,3 in long format and save
////////////////////////////////////////////////////////////////////////////////

keep comp d hr1 lb1 ub1 pval1 hr2 lb2 ub2 pval2
rename d n_pts
rename comp Xvar
reshape long hr lb ub pval, i(Xvar)
rename _j model
rename pval pvalue
rename lb ll
rename ub ul
order model Xvar n_pts hr ll ul
drop if Xvar == "ATCMR"
drop if Xvar == "ABMR"
* list, sepby(Xvar)
gen str21 Author = "Park - 2019"
save park_long, replace
restore

////////////////////////////////////////////////////////////////////////////////
**#  End preparing Park's 2019 Dataset
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
**# Prepare & save Crude and Adj Est & 95% CI from Alachkar/Bagnasco 2022 
////////////////////////////////////////////////////////////////////////////////

clear
cd "C:\Documenti\Delsante\Bagnasco IgA"
use NEW58_IgAN_merged
stset TimefromBiopsytolastFUYr, fail(Graftfailure1Nofailure0 == 1)

discard
postfile result model str21 Xvar n_pts hr ll ul pvalue using fileondisk0, replace
foreach var in  M E S C T {
	preserve
	di _newline(3) in wh "--- CRUDE Analysis of `var'--"
	qui stcox  `var' 
	matrix cstats = r(table)
	local n_pts = e(N) 
	qui test `var' = 0
	local pvalue = r(p)
	post result (`=1') ("`var'") (`n_pts') (`=el(cstats,1,1)')  (`=el(cstats,5,1)') (`=el(cstats,6,1)')  (`=`pvalue'')
	
	di _newline(3) in wh "--- ADJUTED Analysis of `var'--"
	qui stcox  `var' Ageatbiopsy GenderF1M0 YearsfollowuptimefromTXto LN_eGFR_at_biopsy RejectionposttransplantYes1
	matrix astats = r(table)
	local n_pts = e(N) 
	qui test `var' = 0
	local pvalue = r(p)
	post result (`=2') ("`var'") (`n_pts') (`=el(astats,1,1)')  (`=el(astats,5,1)') (`=el(astats,6,1)')  (`=`pvalue'')
	restore
	 }
postclose result



discard
postfile result model str21 Xvar n_pts hr ll ul pvalue using fileondisk1, replace
foreach var in  T1T2 {
	preserve
	di _newline(3) in wh "--- CRUDE Analysis of `var'--"
	qui stcox  i.`var' 
	matrix cstats = r(table)
	local n_pts = e(N) 
	qui test 1.`var' = 0
	local pvalue1 = r(p)
	qui test 2.`var' = 0
	local pvalue2 = r(p)
	post result (`=1') ("T1") (`n_pts') (`=el(cstats,1,2)')  (`=el(cstats,5,2)') (`=el(cstats,6,2)')  (`=`pvalue1'')
	post result (`=1') ("T2") (`n_pts') (`=el(cstats,1,3)')  (`=el(cstats,5,3)') (`=el(cstats,6,3)')  (`=`pvalue2'')
	
	di _newline(3) in wh "--- ADJUTED Analysis of `var'--"
	qui stcox  i.`var' Ageatbiopsy GenderF1M0 YearsfollowuptimefromTXto LN_eGFR_at_biopsy RejectionposttransplantYes1
	matrix astats = r(table)
	local n_pts = e(N) 
	qui test 1.`var' = 0
	local pvalue1 = r(p)
	qui test 2.`var' = 0
	local pvalue2 = r(p)
	post result (`=2') ("T1") (`n_pts') (`=el(astats,1,2)')  (`=el(astats,5,2)') (`=el(astats,6,2)')  (`=`pvalue1'')
	post result (`=2') ("T2") (`n_pts') (`=el(astats,1,3)')  (`=el(astats,5,3)') (`=el(astats,6,3)')  (`=`pvalue2'')
	restore
	 }
postclose result


preserve
clear
use fileondisk0
gen Author = "Bagnasco/Alachkar - 2022"
save fileondisk0, replace
restore

preserve
clear
use fileondisk1
gen Author = "Bagnasco/Alachkar - 2022"
save fileondisk1, replace
restore


preserve
use fileondisk0, clear
append using fileondisk1
append using park_long
* list
foreach var of varlist hr ll ul  {
		gen log`var' = log(`var')
		 }
save bagnasco_park_appended, replace
restore

///////////////////////////////////////////////////////////////////////////////
**# Start Meta-Analysis of the Crude and Adjusted Estimates and export results
//////////////////////////////////////////////////////////////////////////////

putdocx clear 
putdocx begin
putdocx paragraph, style(Title)
putdocx text ("Meta-Analysis with Park's 2019 Study")

preserve
use bagnasco_park_appended, clear
label define model 1 "Crude" 2 "Adjusted"
label values model model
foreach name in M E S C T1 T2  {
	foreach num of numlist 1 2  {
	di _newline(3) in ye "----------------------------------> Component `name' - `: label (model) `num'' Model"
	meta set loghr logll logul if Xvar == "`name'" & model == `num', civartolerance(1e-1) studylabel(Author)
	meta summarize, eform
	meta forestplot, eform("HR") title("`: label (model) `num'' Hazard Ratio (HR)") name(meta_`name'_`num', replace)
	
	putdocx paragraph, style(Heading1)
    putdocx text ("Component `name'")
	
	graph export meta_`name'_`num'.png, replace
	graph export meta_`name'_`num'.tif, replace
	graph export meta_`name'_`num'.pdf, replace
	local hazard_ratio = exp(r(theta)) 
	
	putdocx paragraph, halign(center)
    putdocx image meta_`name'_`num'.png
	 }
	}
restore


preserve
use bagnasco_park_appended, clear
label define model 1 "Crude" 2 "Adjusted"
label values model model
replace Xvar = "C2" if Author == "Bagnasco/Alachkar - 2022" & Xvar == "C"
foreach name in C2  {
	foreach num of numlist 1 2  {
	di _newline(3) in ye "----------------------------------> Component `name' - `: label (model) `num'' Model"
	meta set loghr logll logul if Xvar == "`name'" & model == `num', civartolerance(1e-1) studylabel(Author)
	meta summarize, eform
	meta forestplot, eform("HR") title("`: label (model) `num'' Hazard Ratio (HR)") name(meta_`name'_`num', replace)
	
	putdocx paragraph, style(Heading1)
    putdocx text ("Component `name'")
	
	graph export meta_`name'_`num'.png, replace
	graph export meta_`name'_`num'.tif, replace
	graph export meta_`name'_`num'.pdf, replace
	local hazard_ratio = exp(r(theta)) 
	
	putdocx paragraph, halign(center)
    putdocx image meta_`name'_`num'.png
	
	 }
	}
restore


putdocx paragraph, halign(both)
putdocx text ("Legend of the plots."), bold 
putdocx text (" Random-effects meta-analysis with the findings from Park. ")
putdocx text ("Adjusted model: adjusted for age at biopsy, sex, graft age, eGFR at biopsy, and previous rejection.")
putdocx text (" The values of τ")
putdocx text ("2"), script(super)
putdocx text (", I")
putdocx text ("2"), script(super)
putdocx text (", and H")
putdocx text ("2"), script(super)
putdocx text (" are measures of difference between the HRs (heterogeneity): I")
putdocx text ("2"), script(super) 
putdocx text (" estimates the proportion of variation between the HRs due to heterogeneity (i.e. actual different findings between the two studies) relative to the pure sampling variation. I")
putdocx text ("2"), script(super)
putdocx text (" > 50% indicates substantial heterogeneity; a value of H")
putdocx text ("2"), script(super)
putdocx text (" = 1 indicates perfect homogeneity among the two studies; tau")
putdocx text ("2"), script(super) 
putdocx text (" is the between-study variance.")
putdocx text (" The statistics Q(1) that θ")
putdocx text ("i"), script(sub)
putdocx text (" = θ")
putdocx text ("j"), script(sub)
putdocx text (" is the test that the two HRs are different: a non-significant P value is interpreted as indicating that two HRs are identical, the observed differences being related only to sampling variation.")
putdocx text (" The statistics z that θ = 0 is the test the pooled estimate of two HRs is statistically significant: a significant P value is interpreted as indicating that the histological component is significantly associated with graft failure.")
putdocx text (" The Weight (%) is the proportional contribution of each study to the pooled estimate of the HR (green diamond).")
putdocx text (" The blue square symbol represents the study HR estimate (the dimension of the square symbol is proportional to the Weight of the study); the blue horizontal line represents the associated 95% confidence interval. If the 95% confidence interval does not include the null value (HR=1.0), then the HR is statistically significant (P<0.05).")
putdocx text (" The diamond symbol represents the pooled HR estimate, the center representing the HR, the width the 95% confidence interval")

putdocx save `"meta_analysis_with_Park_`c(current_date)'"', replace
