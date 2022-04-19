capture program drop analyze_TS_completion
capture program drop analyze_TT_completion
capture program drop analyze_TS_burden
capture program drop analyze_TT_burden
capture program drop tabulation_to_matrix

// using round5 data extraction
// 7 clinics in total

program define tabulation_to_matrix, rclass
	local bi_outcome `1'
	local cate_var `2'
	local cond_var `3'
	local condition `4'
	
	dis "***** `condition' ******"
	
	tab `cate_var'
	scalar col_num = r(r)
	
	// save results to matrix
	matrix onerow = J(1,col_num,0)
	matrix rownames onerow = "`condition'"
	
	/// completion rate
	tab `bi_outcome' `cate_var' if `cond_var' == "`condition'" , col matcell(freq) matrow(names)
		
	if r(N) != . {
		if rowsof(freq) == 1 {
			forvalues j = 1/`=col_num' {
					matrix onerow[1,`j'] = names[1,1]
			}
		}
		else {
			forvalues i = 1/1 {
				forvalues j = 1/`=col_num' {
					matrix onerow[`i',`j']= freq[2,`j']/ (freq[1,`j'] + freq[2,`j'])
				}
			}
		}
	}
	
	matrix list onerow
	
	return matrix outcome_rate = onerow
end


// rooming BPA or tobacco screening alert
// effectiveness (encounter level, patient level)
program define analyze_TS_completion
	local enc_type `1'
	// all  vs. selected
	
	use input/alert_enc.dta, clear

	
	// keep TS alerts
	tab alert_desc
	dis "=== physician or resident handled alert ==="
	tab en_alert_desc if prov_type=="Physician" | prov_type=="Resident" 
	drop if !regexm(alert_desc, "PSH.*ROOMING")
	
	// check/remove duplicate alert instances
	// i.e., an alert instance that had two actions: activity link & acknowledge/override
	unique alt_csn_id
	if r(N) > r(unique) {
		dis "=== keep only one action for each alert instance ==="
		log off
		sort alt_csn_id
		quietly by alt_csn_id:  gen dup_alt = cond(_N==1,0,_n)
		sort alt_csn_id
		list alt_csn_id dup_alt trigger_action action_taken prov_type if dup_alt != 0
		list alt_csn_id dup_alt trigger_action en_action_taken if dup_alt != 0 & en_action_taken != 2
		log on
		unique alt_csn_id if dup_alt != 0
		by alt_csn_id: drop if dup_alt != 0 & en_action_taken != 2	
		drop dup_alt
	}
	unique alt_csn_id
	
	
	// distribution of trigger actions vs. action taken
	tab trigger_action, missing
	
	dis "=== action for General BPA section ==="
	tab en_action_taken if trigger_action == "General BPA section", missing
	
	dis "=== action for Open Patient Chart ==="
	tab en_action_taken if trigger_action == "Open Patient Chart", missing
	
	gen time_cate=.
	replace time_cate = 0 if alt2chkin_hr < 0.5
	replace time_cate = 1 if alt2chkin_hr >= 0.5 & alt2chkin_hr < 1
	replace time_cate = 2 if alt2chkin_hr >= 1 & alt2chkin_hr < 2
	replace time_cate = 3 if alt2chkin_hr >= 2 & alt2chkin_hr != .
	label define time_cate_l 0 "<0.5hr" 1 "0.5-1hr" 2 "1-2hr" 3 "2hr later"
	label value time_cate time_cate_l
	
	gen screened_i = 0
	replace screened_i = 1 if override_reason == "Yes" & en_action_taken == 2 & trigger_action != "General BPA section"
		
	gen screened_ni = 0
	replace screened_ni = 1 if override_reason == "Yes" & en_action_taken == 2 & trigger_action == "General BPA section"
	
	gen screened = screened_i + screened_ni
	
	// encoding encounter type
	encode enc_type, gen(en_enc_type)
	tab enc_type en_enc_type, nolabel
	// 6: initial consultation 9: office visit  16: return visit
	// 8: nurse only  11: Post-Op  12: Procedure visit
	// 17: Telemedicine 19: treatment 

	// grouping of encounter types
	gen sel_enc_grp1 = 0
	replace sel_enc_grp1 = 1 if en_enc_type == 6 | en_enc_type == 9 | en_enc_type == 16 ///
	| en_enc_type == 8 | en_enc_type == 11 | en_enc_type == 17 | en_enc_type == 19 ///
	| en_enc_type == 12
	
	save input/alertTS_enc.dta, replace
	
	dis "=== encounter level alert effectiveness ==="
	use input/alertTS_enc.dta, clear
	sort encounter_num 
	collapse (count) num_alt=en_trigger_action (first) user_id department_name enc_type ///
	contact_num dept_abbreviation sel_enc_grp* ///
	en_action_taken period period_variable patient_num ///
	(max) screened screened_i screened_ni, by(encounter_num)
	label value period peri_l
	
	merge m:1 dept_abbreviation using input/clinic_factors.dta
	drop _merge
	
	tab clinic period
	
		
	// shorten department name
	gen department_name_orig=department_name
	replace department_name = regexr(department_name, "(ONCOLOGY)", "")
	replace department_name = regexr(department_name, "(TOLOGY)", "")
	replace department_name = regexr(department_name, "(TION)", "")
	replace department_name = regexr(department_name, "(VIVORSHIP)", "")

	// alert completion by clinics
	dis "=== TS effectiveness (all encounter types) === "
	tab department_name screened , row
	
	dis "=== TS effectiveness (selected encounter types) === "
	tab department_name screened if sel_enc_grp1 == 1, row
	
	
	if "`enc_type'" == "allEnc" {
		dis "=== keep all encounter types ==="
	}
	else if "`enc_type'" == "selEnc" {
		dis "=== keep selected encounter types ==="
		keep if sel_enc_grp1 == 1
		tab enc_type
	}
	
	save input/alertTS_enc_`enc_type'.dta, replace

	// alert completion over implementation time
	use input/alertTS_enc_`enc_type'.dta, clear
	tab screened period, col
	
	tab period
	scalar period_num = r(r)
	matrix ts_enc = J(1,period_num,0)
	levelsof period, local(period_ls)
	matrix colnames ts_enc = `period_ls' 
	matrix rownames ts_enc = clinic
	
	levelsof clinic, local(clinic_ls) 
	foreach c of local clinic_ls {
		/// number of encounters
		tabulation_to_matrix screened period clinic "`c'"
		matrix ts_enc = ts_enc \ r(outcome_rate)
	}
	
	putexcel set "TS_results.xlsx", sheet("TS_eff_`enc_type'") modify
	putexcel A1=matrix(ts_enc), names

	sort period
	collapse (count) num_enc = screened (sum) num_screened = screened, by(period)
	gen rate = num_screened / num_enc
	order period rate num_screened num_enc
	list
	export excel using "TS_results.xlsx", sheet("TS_eff_impltime_`enc_type'") sheetreplace firstrow(variables)
	
	dis "=== patient level (12-month) alert effectiveness ==="
	use input/alertTS_enc_`enc_type'.dta, clear
	foreach var in user_id clinic enc_type period encounter_num {
		local newname = substr("`var'", 1, 8)
		egen tag_`newname' = tag(`var' patient_num)
		egen dist_`newname' = total(tag_`newname'), by(patient_num)
	}
	
	sort patient_num encounter_num
	collapse (count) num_enc=encounter_num ///
	(first) dist_user_id dist_clinic dist_enc_type dist_period dept_abbreviation clinic ///
	(max) screened screened_i screened_ni dist_encount, by(patient_num)
	
	foreach v of var * { 
		label var `v' "" 
	} 

	tab num_enc
	tab dist_period
	tab dist_encoun
	
	// alert completion by clinics
	tab clinic screened if dist_clinic == 1, row
	
	encode clinic, gen(en_clinic)
	
	// alert completion by patient factors
	label define screenedl 0 "not screened" 1 "screened"
	label value screened screenedl 
	
	merge 1:1 patient_num using input/patient_factors.dta
	
	gen sex_race = race_cat3 + (en_sex-1)*3
	label define sex_race_l  0 "F-Other" 1 "F-African American" 2 "F-White"   ///
	3 "M-Other" 4 "M-African American" 5 "M-White"
	label value sex_race sex_race_l

	
	dis "=== all patients ==="
	foreach iv in age_cat3 sex race_cat3 sex_race {
		tab `iv' screened , row chi exact
	}
	
	logistic screened i.age_cat3 en_sex ib(2).race_cat3 i.en_clinic
	logistic screened i.age_cat3 ib(5).sex_race i.en_clinic
	

	dis "=== patients with race White or African American ==="
	foreach iv in age_cat3 sex race_cat3 sex_race {
		tab `iv' screened if race_cat3 != 0, row chi exact
	}
	
	logistic screened i.age_cat3 en_sex ib(2).race_cat3 i.en_clinic if race_cat3 != 0
	logistic screened i.age_cat3 ib(5).sex_race i.en_clinic if race_cat3 != 0
	
end


// treatment BPA or tobacco referral alert
// effectiveness (encounter level, patient level)
// using allEnc, as the encounter types having TT alerts all belong to selEnc types  
program define analyze_TT_completion
	local enc_type allEnc
	
	use input/alert_enc.dta, clear

	
	// keep TT alerts
	tab alert_desc
	dis "=== physician or resident handled alert ==="
	tab en_alert_desc if prov_type=="Physician" | prov_type=="Resident" 
	drop if regexm(alert_desc, "PSH.*ROOMING")
	
	// record TT alert response actions
	gen pp = 0
	replace pp = 1 if override_reason == "Postponed"
	
	gen discussed = 0
	replace discussed = 1 if override_reason == "Initiating Protocol"  | override_reason == "Yes"
	
	gen readyquit = 0
	replace readyquit = 1 if override_reason == "Yes"
	
	gen referral = 0
	replace referral = 1 if order_name == "Amb Referral to Tobacco Cessation"
	
	// resolving inconsistency
	replace readyquit = 1 if referral == 1
	replace discussed = 1 if referral == 1
	replace pp = 0 if referral == 1
	
	
	sort alt_csn_id
	collapse (first) encounter_num user_id department_name enc_type prov_type ///
	contact_num dept_abbreviation *trigger_action override_reason ///
	action_taken en_action_taken period period_variable patient_num alt2chkin_hr ///
	(max) pp discussed readyquit referral, by(alt_csn_id)
	label value period peri_l
	
	
	// distribution of trigger actions vs. action taken
	tab trigger_action, missing
	
	dis "=== action for General BPA section ==="
	tab en_action_taken if trigger_action == "General BPA section", missing
	
	dis "=== action for Open Patient Chart ==="
	tab en_action_taken if trigger_action == "Open Patient Chart", missing

	gen referral_i = 0
	replace referral_i = 1 if referral == 1 & trigger_action != "General BPA section"
		
	gen referral_ni = 0
	replace referral_ni = 1 if referral == 1 & trigger_action == "General BPA section"
	
	// encoding encounter type
	encode enc_type, gen(en_enc_type)
	tab enc_type en_enc_type, nolabel
	// 1: initial consult, 2: Office Visit, 3: Post-Op
	// 4: Procedure visit, 5: Return Patient, 6: Treatment
	
	save input/alertTT_enc.dta, replace
	
	dis "=== encounter level alert effectiveness ==="
	use input/alertTT_enc.dta, clear
	sort encounter_num 
	collapse (count) num_alt=en_trigger_action (first) user_id department_name enc_type ///
	prov_type contact_num dept_abbreviation  ///
	en_action_taken period period_variable patient_num ///
	(max) pp discussed readyquit referral*, by(encounter_num)
	label value period peri_l
	
	merge m:1 dept_abbreviation using input/clinic_factors.dta
	drop if _merge != 3
	drop _merge
	
	// resolving encounter-level inconsistency
	replace readyquit = 1 if referral == 1
	replace discussed = 1 if readyquit == 1
	replace pp = 0 if discussed == 1
	
	tab clinic period
	
	// shorten department name
	gen department_name_orig=department_name
	replace department_name = regexr(department_name, "(ONCOLOGY)", "")
	replace department_name = regexr(department_name, "(TOLOGY)", "")
	replace department_name = regexr(department_name, "(TION)", "")
	replace department_name = regexr(department_name, "(VIVORSHIP)", "")

	// referral completion by clinics
	dis "=== TT effectiveness (all clinics) === "
	tab clinic referral , row
		
	if "`enc_type'" == "allEnc" {
		dis "=== keep all encounter types ==="
	}
	else if "`enc_type'" == "selEnc" {
		dis "=== keep selected encounter types ==="
		// not implement yet
	}
	
	save input/alertTT_enc_`enc_type'.dta, replace

	// alert completion
	use input/alertTT_enc_`enc_type'.dta, clear
	collapse pp discussed readyquit referral
	replace pp = pp * 100
	gen npp = 100 - pp
	replace discussed = discussed * 100
	replace readyquit = readyquit * 100
	replace referral = referral * 100
	list
	drop pp
	order npp discussed readyquit referral
	export excel using "TT_results.xlsx", sheet("TT_var_eff_`enc_type'") sheetreplace firstrow(variables)

		
	// alert completion over implementation time
	use input/alertTT_enc_`enc_type'.dta, clear
	tab referral period, col
	
	tab period
	scalar period_num = r(r)
	matrix tt_enc = J(1,period_num,0)
	levelsof period, local(period_ls)
	matrix colnames tt_enc = `period_ls' 
	matrix rownames tt_enc = clinic
	
	levelsof clinic, local(clinic_ls) 
	foreach c of local clinic_ls {
		/// number of encounters
		tabulation_to_matrix referral period clinic "`c'"
		matrix tt_enc = tt_enc \ r(outcome_rate)
	}
	
	putexcel set "TT_results.xlsx", sheet("TT_eff_`enc_type'") modify
	putexcel A1=matrix(tt_enc), names

	sort period
	collapse (count) num_enc = referral (sum) num_referral = referral, by(period)
	gen rate = num_referral / num_enc
	order period rate num_referral num_enc
	list
	export excel using "TT_results.xlsx", sheet("TT_eff_impltime_`enc_type'") sheetreplace firstrow(variables)
	
	
	dis "=== patient level (12-month) alert effectiveness ==="
	use input/alertTT_enc_`enc_type'.dta, clear
	foreach var in user_id clinic enc_type period encounter_num {
		local newname = substr("`var'", 1, 8)
		egen tag_`newname' = tag(`var' patient_num)
		egen dist_`newname' = total(tag_`newname'), by(patient_num)
	}
	
	sort patient_num encounter_num
	collapse (count) num_enc=encounter_num ///
	(first) dist_user_id dist_clinic dist_enc_type dist_period dept_abbreviation clinic ///
	(max) pp discussed readyquit referral* dist_encount, by(patient_num)
		
	foreach v of var * { 
		label var `v' "" 
	} 

	tab num_enc
	tab dist_period
	tab dist_encoun
	
	// alert completion by clinics
	tab clinic referral if dist_clinic == 1, row
	
	
	// alert completion by patient factors
	label define referrall 0 "not referred" 1 "referred"
	label value referral referrall 
	
	merge 1:1 patient_num using input/patient_factors.dta
	gen sex_race = race_cat3 + (en_sex-1)*3
	label define sex_race_l  0 "F-Other" 1 "F-African American" 2 "F-White"   ///
	3 "M-Other" 4 "M-African American" 5 "M-White"
	label value sex_race sex_race_l
	
	dis "=== all patients ==="
	foreach iv in age_cat3 sex race_cat3 sex_race {
		tab `iv' referral, row chi exact
	}
	
	foreach iv in age_cat3 sex race_cat3 sex_race {
		tab `iv' discussed, row chi exact
	}
	
	dis "=== patients with race White or African American ==="
	foreach iv in age_cat3 sex race_cat3 sex_race {
		tab `iv' referral if race_cat3 != 0, row chi exact
	}
	
	foreach iv in age_cat3 sex race_cat3 sex_race {
		tab `iv' discussed if race_cat3 != 0, row chi exact
	}
end


// alert burden, interruptive alerts only
// run this code after running analyze_TS_completion
program define analyze_TS_burden
	local enc_type `1'
	// all vs. selected encounters
	
	local enc_type allEnc
	*** alert burden metrics 1 : number to fire ***
	use input/alertTS_enc.dta, clear
	drop if trigger_action == "General BPA section"
	merge m:1 dept_abbreviation using input/clinic_factors.dta
	drop _merge
	
	replace screened = screened_i
	gen time4alert = ack_tkn_inst_s2i - alt_action_inst_s2i
	gen time4alert_pp = time4alert
	replace time4alert_pp = 0 if screened == 1
	gen time4alert_np = time4alert - time4alert_pp
	dis "=== descriptive statistics of alert response time ==="
	dis "== per alert instance (all alert instances) =="
	tabstat time4alert_np time4alert_pp time4alert, stat(mean, sd, p50)
	
	dis "== per alert instance (instances not postponed) =="
	tabstat time4alert_np if screened == 1, stat(mean, sd, p50, min, max)
	
	dis "== per alert instance (instances postponed) =="
	tabstat time4alert_pp if screened == 0, stat(mean, sd, p50, min, max)
	
	// shorten department name
	gen department_name_orig=department_name
	replace department_name = regexr(department_name, "(ONCOLOGY)", "")
	replace department_name = regexr(department_name, "(TOLOGY)", "")
	replace department_name = regexr(department_name, "(TION)", "")
	replace department_name = regexr(department_name, "(VIVORSHIP)", "")
	
	
	tab clinic period
		
	*** alert burden metrics 1: alert level (number to fire); 2, 3, 4 : encounter level ***
	sort encounter_num 
	collapse (count) num_alt=en_trigger_action (first) user_id department_name_orig department_name ///
	clinic enc_type contact_num dept_abbreviation sel_enc_grp* ///
	en_action_taken period period_variable patient_num ///
	(sum) time4alert_pp = time4alert_pp time4alert = time4alert screened screened_i ///
	, by(encounter_num)
	label value period peri_l
	
	
	if "`enc_type'" == "allEnc" {
		dis "=== keep all encounter types ==="
	}
	else if "`enc_type'" == "selEnc" {
		dis "=== keep selected encounter types ==="
		keep if sel_enc_grp == 1
		tab enc_type
	}
	
	save input/alertTS_enc_burden_`enc_type'.dta, replace

	// alert burden by clinics
	use input/alertTS_enc_burden_`enc_type'.dta, clear
	sort clinic
	gen time4alert_np =  time4alert - time4alert_pp
	foreach tvar in time4alert_np time4alert_pp time4alert {
		dis
		dis "sum `tvar'"
		sum `tvar'
		
		dis
		dis "sum `tvar' if `tvar' != 0"
		sum `tvar' if `tvar' != 0
	}
	
	list screened screened_i time4alert_np time4alert_pp if time4alert == 0
	
	sum num_alt
	scalar num_alt_fired = r(sum)
	sum screened
	scalar num_alt_completed = r(sum)
	dis "number to fire", num_alt_fired/num_alt_completed
	
	
	collapse (count) num_enc = encounter_num ///
	(mean) time4alert_pp_mean = time4alert_pp time4alert_mean = time4alert ///
	num_alt_mean = num_alt num_rs_mean = screened ///
	(sd) time4alert_pp_sd = time4alert_pp time4alert_sd = time4alert ///
	num_alt_sd = num_alt num_rs_sd = screened ///
	(sum) time4alert_pp = time4alert_pp time4alert = time4alert ///
	num_alt = num_alt num_rs = screened ///
	, by(clinic)
	
	gen num2fire = num_alt / num_rs
	gen pp_time_rate = time4alert_pp / time4alert
	gen time4alert_np_mean = time4alert_mean - time4alert_pp_mean
	
	list
	keep clinic num2fire pp_time_rate num_alt_mean num_alt_sd time4alert_pp_mean time4alert_np_mean time4alert_mean
	order clinic num2fire pp_time_rate num_alt_mean num_alt_sd time4alert_pp_mean time4alert_np_mean time4alert_mean
	export excel using "TS_results.xlsx", sheet("TS_burden_`enc_type'") sheetreplace firstrow(variables)
	
	
end

//using allEnc, as the encounter types having TT alerts all belong to selEnc types  
program define analyze_TT_burden
	local enc_type allEnc
	
	use input/alert_enc.dta, clear
	
	// keep TT alerts
	tab alert_desc
	dis "=== physician or resident handled alert ==="
	tab en_alert_desc if prov_type=="Physician" | prov_type=="Resident" 
	drop if regexm(alert_desc, "PSH.*ROOMING")
	drop if trigger_action == "General BPA section"
	
	// record TT alert response actions
	gen pp = 0
	replace pp = 1 if override_reason == "Postponed"
	
	gen discussed = 0
	replace discussed = 1 if override_reason == "Initiating Protocol"  | override_reason == "Yes"
	
	gen readyquit = 0
	replace readyquit = 1 if override_reason == "Yes"
	
	gen referral = 0
	replace referral = 1 if order_name == "Amb Referral to Tobacco Cessation"
	
	
	// resolving inconsistency
	replace readyquit = 1 if referral == 1
	replace discussed = 1 if referral == 1
	replace pp = 0 if referral == 1

	sort alt_csn_id 
	collapse (first) encounter_num user_id department_name enc_type prov_type ///
	contact_num dept_abbreviation *trigger_action override_reason ///
	action_taken en_action_taken period period_variable patient_num alt2chkin_hr ///
	(min) alt_action_inst_s2i ///
	(max) pp discussed readyquit referral ack_tkn_inst_s2i, by(alt_csn_id)
	label value period peri_l
	
	gen time_cate=.
	replace time_cate = 0 if alt2chkin_hr < 0.5
	replace time_cate = 1 if alt2chkin_hr >= 0.5 & alt2chkin_hr < 1
	replace time_cate = 2 if alt2chkin_hr >= 1 & alt2chkin_hr < 2
	replace time_cate = 3 if alt2chkin_hr >= 2 & alt2chkin_hr != .
	label define time_cate_l 0 "<0.5hr" 1 "0.5-1hr" 2 "1-2hr" 3 "2hr later"
	label value time_cate time_cate_l
	
	// encoding encounter type
	encode enc_type, gen(en_enc_type)
	tab enc_type en_enc_type, nolabel
	// 1: initial consult, 2: Office Visit, 3: Post-Op
	// 4: Procedure visit, 5: Return Patient, 6: Treatment
	
	gen time4alert = ack_tkn_inst_s2i - alt_action_inst_s2i
	gen time4alert_pp = time4alert
	replace time4alert_pp = 0 if pp == 0
	gen time4alert_np = time4alert - time4alert_pp
	dis "=== descriptive statistics of alert response time ==="
	dis "== per alert instance (all alert instances) =="
	tabstat time4alert_np time4alert_pp time4alert, stat(mean, sd, p50)
	
	dis "== per alert instance (instances not postponed) =="
	tabstat time4alert_np if pp == 0, stat(mean, sd, p50, min, max)
	
	dis "== per alert instance (instances postponed) =="
	tabstat time4alert_pp if pp == 1, stat(mean, sd, p50, min, max)
	
	// shorten department name
	gen department_name_orig=department_name
	replace department_name = regexr(department_name, "(ONCOLOGY)", "")
	replace department_name = regexr(department_name, "(TOLOGY)", "")
	replace department_name = regexr(department_name, "(TION)", "")
	replace department_name = regexr(department_name, "(VIVORSHIP)", "")
	save input/alertTT_burden_enc.dta, replace
	
	*** alert burden metrics 1: alert level (number to fire); 2, 3, 4 : encounter level ***	
	use input/alertTT_burden_enc.dta, clear
	gen np = 1 - pp
	sort encounter_num 
	collapse (count) num_alt=en_trigger_action (first) user_id department_name enc_type ///
	prov_type contact_num dept_abbreviation  ///
	en_action_taken period period_variable patient_num ///
	(sum) time4alert_pp = time4alert_pp time4alert = time4alert ///
	pp np discussed readyquit referral, by(encounter_num)
	label value period peri_l
	
	merge m:1 dept_abbreviation using input/clinic_factors.dta
	drop if _merge != 3
	drop _merge
	
	if "`enc_type'" == "allEnc" {
		dis "=== keep all encounter types ==="
	}
	else if "`enc_type'" == "selEnc" {
		dis "=== keep selected encounter types ==="
		keep if sel_enc_grp1 == 1
		tab enc_type
	}
	
	save input/alertTT_enc_burden_`enc_type'.dta, replace

	// alert burden by clinics
	use input/alertTT_enc_burden_`enc_type'.dta, clear
	sort clinic
	gen time4alert_np =  time4alert - time4alert_pp
	foreach tvar in time4alert_np time4alert_pp time4alert {
		dis
		dis "sum `tvar'"
		sum `tvar'
		
		dis
		dis "sum `tvar' if `tvar' != 0"
		sum `tvar' if `tvar' != 0
	}
	
	list referral time4alert_np time4alert_pp if time4alert == 0
	
	
	collapse (count) num_enc = encounter_num ///
	(mean) time4alert_pp_mean = time4alert_pp time4alert_mean = time4alert ///
	num_alt_mean = num_alt num_rs_mean = np ///
	(sd) time4alert_pp_sd = time4alert_pp time4alert_sd = time4alert ///
	num_alt_sd = num_alt num_pp_sd = np ///
	(sum) time4alert_pp = time4alert_pp time4alert = time4alert ///
	num_alt = num_alt num_pp = pp ///
	, by(clinic)
	
	gen num2fire = num_alt / (num_alt - num_pp)
	gen pp_time_rate = time4alert_pp / time4alert
	gen time4alert_np_mean = time4alert_mean - time4alert_pp_mean
	list
	keep clinic num2fire pp_time_rate num_alt_mean num_alt_sd time4alert_pp_mean time4alert_np_mean time4alert_mean
	order clinic num2fire pp_time_rate num_alt_mean num_alt_sd time4alert_pp_mean time4alert_np_mean time4alert_mean
	export excel using "TT_results.xlsx", sheet("TT_burden_`enc_type'") sheetreplace firstrow(variables)	
end



//*
foreach enc_type in allEnc selEnc { // selEnc allEnc 
	log using "logfiles/analyze_TS_completion_`enc_type'.log", replace
	analyze_TS_completion `enc_type'
	log close
}
*/

//*
log using "logfiles/analyze_TT_completion_allEnc.log", replace
analyze_TT_completion
log close
*/

//*
foreach enc_type in selEnc allEnc { // selEnc allEnc 
	log using "logfiles/analyze_TS_burden_`enc_type'.log", replace
	analyze_TS_burden `enc_type'
	log close
}
*/


//*
log using "logfiles/analyze_TT_burden_allEnc.log", replace
analyze_TT_burden
log close
*/


