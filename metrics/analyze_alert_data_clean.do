***************************************
* MIT License
* Copyright (c) 2023 Jinying Chen
*  
* author: Jinying Chen, iDAPT Cancer Control Center
* date: 2023-2-20
* ver: 2.0
* 
* This code was written to support data analysis for the MIER project and 
* the 2023 paper published in JMIR Medical Informatics.
* The code is for research use only, and is provided as it is.
* 

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
	
	
	// load encounter level EHR documentation for TS results
	use input/screening_flowsheet.dta, clear
	keep encounter_num flo* sc_result
	save input/screening_flowsheet_results.dta, replace
	
	use input/alert_enc.dta, clear
	
	// keep TS alerts
	tab alert_desc
	dis "=== physician or resident handled alert ==="
	tab en_alert_desc if prov_type=="Physician" | prov_type=="Resident" 
	drop if !regexm(alert_desc, "PSH.*ROOMING")
	
	// handle multiple actions
	// an alert instance may have 1-3 actions: activity link, acknowledge/override, Send In Basket Message
	gen activityLink = 0
	gen acknowledge = 0
	gen sendBasketMsg = 0
	replace activityLink = 1 if regexm(action_taken, "Activity") 
	replace acknowledge = 1 if regexm(action_taken, "Acknowledge")
	replace sendBasketMsg = 1 if regexm(action_taken, "Basket")
	sort encounter_num alt_id alt_csn_id 
	collapse (firstnm) trigger_action override_reason user_id department_name ///
	contact_num dept_abbreviation enc_type ///
	period period_variable patient_num encounter_num alt2chkin_hr ///
	alt_action_inst_s2i ack_tkn_inst_s2i /// 
	(max) activityLink acknowledge sendBasketMsg, by(alt_csn_id)
	
	
	// distribution of trigger actions vs. action taken
	tab trigger_action, missing
	
	gen ignore_alt = 0
	replace ignore_alt = 1 if activityLink + acknowledge + sendBasketMsg == 0
	label define ignore_alt_l 0 "not ignored" 1 "ignored"
	label value ignore_alt ignore_alt_l
	
	dis "=== action for General BPA section ==="
	tab ignore_alt if trigger_action == "General BPA section", missing
	
	dis "=== action for Open Patient Chart ==="
	tab ignore_alt if trigger_action == "Open Patient Chart", missing
	
	gen time_cate=.
	replace time_cate = 0 if alt2chkin_hr < 0.5
	replace time_cate = 1 if alt2chkin_hr >= 0.5 & alt2chkin_hr < 1
	replace time_cate = 2 if alt2chkin_hr >= 1 & alt2chkin_hr < 2
	replace time_cate = 3 if alt2chkin_hr >= 2 & alt2chkin_hr != .
	label define time_cate_l 0 "<0.5hr" 1 "0.5-1hr" 2 "1-2hr" 3 "2hr later"
	label value time_cate time_cate_l
	
	gen screened_ack_i = 0
	replace screened_ack_i = 1 if override_reason == "Yes" & trigger_action != "General BPA section"
		
	gen screened_ack_ni = 0
	replace screened_ack_ni = 1 if override_reason == "Yes" & trigger_action == "General BPA section"
	
	gen screened_ack = screened_ack_i + screened_ack_ni
	
	
	gen screened_act_i = 0
	replace screened_act_i = 1 if override_reason == "Yes" & activityLink + sendBasketMsg > 0 & trigger_action != "General BPA section"
		
	gen screened_act_ni = 0
	replace screened_act_ni = 1 if override_reason == "Yes" & activityLink + sendBasketMsg > 0 & trigger_action == "General BPA section"
	
	gen screened_act = screened_act_i + screened_act_ni
	
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
	collapse (count) num_alt = alt_csn_id (first) user_id department_name enc_type ///
	contact_num dept_abbreviation sel_enc_grp* ///
	period period_variable patient_num ///
	(max) screened*, by(encounter_num)
	label value period peri_l
	
	merge m:1 dept_abbreviation using input/clinic_factors.dta
	drop _merge
	
	tab clinic period
	
	merge 1:1 encounter_num using input/screening_flowsheet_results.dta
	drop if _merge == 2
	rename _merge merge_flowsheet
	replace flo_s12_cmpl = 0 if flo_s12_cmpl == .
	replace flo_part_cmpl = 0 if flo_part_cmpl == .
	
	// shorten department name
	gen department_name_orig=department_name
	replace department_name = regexr(department_name, "(ONCOLOGY)", "")
	replace department_name = regexr(department_name, "(TOLOGY)", "")
	replace department_name = regexr(department_name, "(TION)", "")
	replace department_name = regexr(department_name, "(VIVORSHIP)", "")

	// alert completion acknowledged by clinics
	dis "=== TS completion acknowledged (all encounter types) === "
	tab department_name screened_ack , row
	
	dis "=== TS completion acknowledged (selected encounter types) === "
	tab department_name screened_ack if sel_enc_grp1 == 1, row
	
	dis "=== TS acted upon (all encounter types) === "
	tab department_name screened_act , row
	
	dis "=== TS acted upon (selected encounter types) === "
	tab department_name screened_act if sel_enc_grp1 == 1, row
	
	
	dis "=== TS complete flowsheet (all encounter types) === "
	tab department_name flo_part_cmpl, row
	
	dis "=== TS complete flowsheet (selected encounter types) === "
	tab department_name flo_part_cmpl if sel_enc_grp1 == 1, row
	
	
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
	foreach screened_var in screened_ack  flo_part_cmpl {
		use input/alertTS_enc_`enc_type'.dta, clear
		gen screened = `screened_var'
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
	
		putexcel set "TS_results_`screened_var'.xlsx", sheet("TS_eff_`enc_type'") modify
		putexcel A1=matrix(ts_enc), names

		sort period
		collapse (count) num_enc = screened (sum) num_screened = screened, by(period)
		gen rate = num_screened / num_enc
			order period rate num_screened num_enc
		list
		export excel using "TS_results_`screened_var'.xlsx", sheet("TS_eff_impltime_`enc_type'") sheetreplace firstrow(variables)
	}
	
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
	(max) screened_* dist_encount, by(patient_num)
	
	foreach v of var * { 
		label var `v' "" 
	} 

	tab num_enc
	tab dist_period
	tab dist_encoun
	
	// alert completion by clinics
	tab clinic screened_ack if dist_clinic == 1, row
	
	encode clinic, gen(en_clinic)
	
	// alert completion by patient factors
	label define screenedl 0 "not screened" 1 "screened"
	label value screened_ack screenedl 
	
	merge 1:1 patient_num using input/patient_factors.dta
	
	gen sex_race = race_cat3 + (en_sex-1)*3
	label define sex_race_l  0 "F-Other" 1 "F-African American" 2 "F-White"   ///
	3 "M-Other" 4 "M-African American" 5 "M-White"
	label value sex_race sex_race_l

	
	dis "=== all patients ==="
	foreach iv in age_cat3 sex race_cat3 sex_race {
		tab `iv' screened_ack , row chi exact
	}
	
	logistic screened_ack i.age_cat3 en_sex ib(2).race_cat3 i.en_clinic
	logistic screened_ack i.age_cat3 ib(5).sex_race i.en_clinic
	

	dis "=== patients with race White or African American ==="
	foreach iv in age_cat3 sex race_cat3 sex_race {
		tab `iv' screened_ack if race_cat3 != 0, row chi exact
	}
	
	logistic screened_ack i.age_cat3 en_sex ib(2).race_cat3 i.en_clinic if race_cat3 != 0
	logistic screened_ack i.age_cat3 ib(5).sex_race i.en_clinic if race_cat3 != 0
	
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
	
	// TT alert response actions
	tab action_taken

	gen openSmartSet = 0
	gen acknowledge = 0
	gen sendBasketMsg = 0
	replace openSmartSet = 1 if regexm(action_taken, "Open SmartSet") 
	replace acknowledge = 1 if regexm(action_taken, "Acknowledge")
	replace sendBasketMsg = 1 if regexm(action_taken, "Basket")
	gen referral = 0
	replace referral = 1 if order_name == "Amb Referral to Tobacco Cessation"
	
	// record TT alert response actions
	// override reasons:
	/*
	a.       Initiate Protocol = Discussed, not ready to quit
	b.       Yes = Discussed, ready to quit
	c.       Clinician Reviewed = Not appropriate (comment required)
	d.       Patient Declines = Not discussed today (defer 24hrs)
	e.       Postponed = Deferred 10 minutes 
	*/
	
	sort alt_csn_id
	collapse (firstnm) encounter_num user_id department_name enc_type prov_type ///
	contact_num dept_abbreviation *trigger_action override_reason ///
	action_taken en_action_taken period period_variable patient_num alt2chkin_hr ///
	alt_action_inst_s2i ack_tkn_inst_s2i /// 
	(max) openSmartSet acknowledge sendBasketMsg referral, by(alt_csn_id)
	label value period peri_l
	
	gen ignore_alt = 0
	replace ignore_alt = 1 if openSmartSet + acknowledge + sendBasketMsg == 0
	label define ignore_alt_l 0 "not ignored" 1 "ignored"
	label value ignore_alt ignore_alt_l
	
		
	//* method used in JMI paper
	gen np = 0
	replace np = 1 if override_reason == "Initiating Protocol"  | override_reason == "Yes" ///
	| override_reason == "Clinician reviewed" |  override_reason == "Patient declines" ///
	| referral > 0
	*/
	
	gen discussed = 0
	replace discussed = 1 if override_reason == "Initiating Protocol"  | override_reason == "Yes"
	
	gen readyquit = 0
	replace readyquit = 1 if override_reason == "Yes"
		
	// resolving inconsistency
	replace readyquit = 1 if referral == 1
	replace discussed = 1 if readyquit == 1
	replace np = 1 if discussed == 1
	
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
	(max) np discussed readyquit referral*, by(encounter_num)
	label value period peri_l
	
	merge m:1 dept_abbreviation using input/clinic_factors.dta
	drop if _merge != 3
	drop _merge
	
	// resolving encounter-level inconsistency
	replace readyquit = 1 if referral == 1
	replace discussed = 1 if readyquit == 1
	replace np = 1 if discussed == 1
	gen pp = 1 - np
	
	tab clinic period
	
	
	// shorten department name
	gen department_name_orig=department_name
	replace department_name = regexr(department_name, "(ONCOLOGY)", "")
	replace department_name = regexr(department_name, "(TOLOGY)", "")
	replace department_name = regexr(department_name, "(TION)", "")
	replace department_name = regexr(department_name, "(VIVORSHIP)", "")

		
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
	foreach iv in pp discussed readyquit referral {
		tab `iv', missing
	}
	
	collapse pp discussed readyquit referral
	replace pp = pp * 100
	gen np = 100 - pp
	replace discussed = discussed * 100
	replace readyquit = readyquit * 100
	replace referral = referral * 100
	
	dis "=== TT effectiveness (overall) ==="
	list
	drop pp
	order np discussed readyquit referral
	export excel using "TT_results.xlsx", sheet("TT_var_eff_`enc_type'") sheetreplace firstrow(variables)

		
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
	(max) np discussed readyquit referral* dist_encount, by(patient_num)
		
	foreach v of var * { 
		label var `v' "" 
	} 
	
	// resolving patient-level inconsistency
	replace readyquit = 1 if referral == 1
	replace discussed = 1 if readyquit == 1
	replace np = 1 if discussed == 1
	gen pp = 1 - np
	
	tab num_enc
	tab dist_period
	tab dist_encoun
	
	
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
	//local enc_type allEnc
	
	use input/alertTS_enc.dta, clear
	
	*** keep interruptive alerts only ***
	drop if trigger_action == "General BPA section"
	merge m:1 dept_abbreviation using input/clinic_factors.dta
	drop _merge
	
	// remove 13 alerts that missed alert handling time point, i.e., ack_tkn_inst_s2i == .
	drop if ignore_alt == 1
	
	gen pp = 0
	replace pp = 1 if regexm(override_reason, "Postponed")
	tab override_reason if pp == 0 & screened_ack_i == 0, missing  // 128 "Does not meet criteria"
	gen np = 1 - pp
	
	gen screened = screened_ack_i
	gen time4alert = ack_tkn_inst_s2i - alt_action_inst_s2i
	
	
	gen time4alert_pp = time4alert
	replace time4alert_pp = 0 if pp == 0
	gen time4alert_np = time4alert - time4alert_pp
	
	// handling outliers
	xtile xt_np = time4alert_np if time4alert_np > 0, nq(20)
	xtile xt_pp = time4alert_pp if time4alert_pp > 0, nq(20)
	sum time4alert_np if xt_np == 19
	scalar time_cap_np1 =  r(max)
	sum time4alert_pp if xt_pp == 19
	scalar time_cap_pp1 =  r(max)
	replace time4alert_np = time_cap_np1  if time4alert_np > time_cap_np1
	replace time4alert_pp = time_cap_pp1  if time4alert_pp > time_cap_pp1
	
	sum time4alert_np if xt_np == 2
	scalar time_cap_np2 =  r(min)
	sum time4alert_pp if xt_pp == 2
	scalar time_cap_pp2 =  r(min)
	replace time4alert_np = time_cap_np2 if xt_np == 1
	replace time4alert_pp = time_cap_pp2 if xt_pp == 1
	replace time4alert = time4alert_np + time4alert_pp
	
	dis "=== total alert handling time ==="
	sum time4alert_np
	dis
	dis "total time spent completing the alert: ", r(sum)/3600
	dis
	
	sum time4alert_pp
	dis
	dis "total time spent postponing the alert: ", r(sum)/3600
	dis
	
	dis "=== alert-level alert handling time ==="
	dis "== per alert instance (acknowledged completion or does not meet criteria) =="
	tabstat time4alert_np if time4alert_np > 0, stat(mean, sd, p50, p25, p75, min, max, n)
	
	dis "== per alert instance (use activity link) =="
	tabstat time4alert_np if screened_act_i == 1, stat(mean, sd, p50, p25, p75, min, max, n)
	
	dis "== per alert instance (instances postponed) =="
	tabstat time4alert_pp if time4alert_pp > 0, stat(mean, sd, p50, p25, p75, min, max, n)
	
	
	// shorten department name
	gen department_name_orig=department_name
	replace department_name = regexr(department_name, "(ONCOLOGY)", "")
	replace department_name = regexr(department_name, "(TOLOGY)", "")
	replace department_name = regexr(department_name, "(TION)", "")
	replace department_name = regexr(department_name, "(VIVORSHIP)", "")
	
	
	tab clinic period
		
	sort encounter_num 
	collapse (count) num_alt=alt_csn_id  (first) user_id department_name_orig department_name ///
	clinic enc_type contact_num dept_abbreviation sel_enc_grp* ///
	period period_variable patient_num ///
	(max) pp = pp np = np ///
	(sum) time4alert_pp = time4alert_pp time4alert = time4alert time4alert_np = time4alert_np screened* ///
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
	
	
	dis "=== encounter-level alert handling time ==="
	foreach tvar in time4alert_np time4alert_pp time4alert {
		dis
		tabstat `tvar' if `tvar' > 0, stat(mean, sd, p50, p25, p75, min, max, n)	
	}
	
	list screened* time4alert_np time4alert_pp if time4alert == 0
	
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
	num_alt = num_alt num_rs = screened num_enc_pp = pp num_enc_np = np ///
	, by(clinic)
	
	gen num2fire = num_alt / num_rs
	gen pp_time_rate = time4alert_pp / time4alert
	gen time4alert_np_mean = time4alert_mean - time4alert_pp_mean
	gen time4alert_np = time4alert - time4alert_pp
	
	list
	keep clinic num2fire pp_time_rate num_alt_mean num_alt_sd ///
	time4alert_pp time4alert_np time4alert num_enc_pp num_enc_np num_enc
	order clinic num2fire pp_time_rate num_alt_mean num_alt_sd ///
	time4alert_pp time4alert_np time4alert num_enc_pp num_enc_np num_enc
	export excel using "TS_burden_results.xlsx", sheet("TS_burden_`enc_type'") sheetreplace firstrow(variables)
	
	
end

//using allEnc, as the encounter types having TT alerts all belong to selEnc types 
//run this code after running analyze_TT_completion
program define analyze_TT_burden
	local enc_type allEnc
	
	use input/alertTT_enc.dta, clear
	
	*** keep interruptive alerts only ***
	drop if trigger_action == "General BPA section"
	
	// remove 7 alerts that were ignored, i.e., ack_tkn_inst_s2i == .
	drop if ignore_alt == 1
	
	tab enc_type en_enc_type, nolabel
	// 1: initial consult, 2: Office Visit, 3: Post-Op
	// 4: Procedure visit, 5: Return Patient, 6: Treatment
	
	gen pp = 1 - np
	
	gen time4alert = ack_tkn_inst_s2i - alt_action_inst_s2i
	gen time4alert_pp = time4alert
	replace time4alert_pp = 0 if pp == 0
	gen time4alert_np = time4alert
	replace time4alert_np = 0 if np == 0
	gen time4alert_ref = time4alert
	replace time4alert_ref = 0 if referral == 0
	gen time4alert_dis = time4alert
	replace time4alert_dis = 0 if discussed == 0
	
	ci means time4alert_np if pp == 1
	ci means time4alert_pp if pp == 0
	
	// handling outliers
	xtile xt_np = time4alert_np if time4alert_np > 0, nq(20)
	xtile xt_pp = time4alert_pp if time4alert_pp > 0, nq(20)
	xtile xt_dis = time4alert_dis if time4alert_dis > 0, nq(20)
	sum time4alert_np if xt_np == 19
	scalar time_cap_np1 =  r(max)
	sum time4alert_pp if xt_pp == 19
	scalar time_cap_pp1 =  r(max)
	sum time4alert_dis if xt_dis == 19
	scalar time_cap_dis1 =  r(max)
	replace time4alert_np = time_cap_np1  if time4alert_np > time_cap_np1
	replace time4alert_pp = time_cap_pp1  if time4alert_pp > time_cap_pp1
	replace time4alert_dis = time_cap_dis1  if time4alert_dis > time_cap_dis1
	
	sum time4alert_np if xt_np == 2
	scalar time_cap_np2 =  r(min)
	sum time4alert_pp if xt_pp == 2
	scalar time_cap_pp2 =  r(min)
	sum time4alert_dis if xt_dis == 2
	scalar time_cap_dis2 =  r(min)
	replace time4alert_np = time_cap_np2 if xt_np == 1
	replace time4alert_pp = time_cap_pp2 if xt_pp == 1
	replace time4alert_dis = time_cap_dis2 if xt_dis == 1
	
	replace time4alert = time4alert_np + time4alert_pp
	
	dis "=== total alert handling time ==="
	sum time4alert_np
	dis
	dis "total time spent completing the alert: ", r(sum)/3600
	dis
	
	sum time4alert_pp
	dis
	dis "total time spent postponing the alert: ", r(sum)/3600
	dis
	
	dis "=== alert-level alert handling time ==="
	dis "== per alert instance (instances postponed) =="
	tabstat time4alert_pp if pp == 1, stat(mean, sd, p50, p25, p75, min, max, n)
	
	dis "== per alert instance (instances not postponed) =="
	tabstat time4alert_np if np == 1, stat(mean, sd, p50, p25, p75, min, max, n)
	
	dis "== per alert instance (instances discussed) =="
	tabstat time4alert_dis if discussed == 1, stat(mean, sd, p50, p25, p75, min, max, n)
	
	dis "== per alert instance (instances resulting a referral) =="
	tabstat time4alert_ref if referral == 1, stat(mean, sd, p50, p25, p75, min, max, n)
	
	// shorten department name
	gen department_name_orig=department_name
	replace department_name = regexr(department_name, "(ONCOLOGY)", "")
	replace department_name = regexr(department_name, "(TOLOGY)", "")
	replace department_name = regexr(department_name, "(TION)", "")
	replace department_name = regexr(department_name, "(VIVORSHIP)", "")
	save input/alertTT_burden_enc.dta, replace
	
	*** alert burden metrics 1: alert level (number to fire); 2, 3, 4 : encounter level ***	
	use input/alertTT_burden_enc.dta, clear
	gen dis_or_ref = 0
	replace dis_or_ref = 1 if discussed == 1 | referral == 1
	
	sort encounter_num 
	collapse (count) num_alt=en_trigger_action (first) user_id department_name enc_type ///
	prov_type contact_num dept_abbreviation  ///
	en_action_taken period period_variable patient_num ///
	(max) enc_pp = pp enc_np = np enc_dis = discussed ///
	(sum) time4alert* pp np discussed readyquit referral dis_or_ref, by(encounter_num)
	label value period peri_l
	
	merge m:1 dept_abbreviation using input/clinic_factors.dta
	drop if _merge == 2
	drop _merge
	
	if "`enc_type'" == "allEnc" {
		dis "=== keep all encounter types ==="
	}
	else if "`enc_type'" == "selEnc" {
		dis "=== keep selected encounter types ==="
		keep if sel_enc_grp3 == 1
		tab enc_type
	}
	
	save input/alertTT_enc_burden_`enc_type'.dta, replace

	// alert burden by clinics
	use input/alertTT_enc_burden_`enc_type'.dta, clear
	sort clinic
	
	dis "=== encounter-level alert handling time ==="
	foreach tvar in time4alert_np time4alert_pp time4alert time4alert_dis {
		dis
		tabstat `tvar' if `tvar' > 0, stat(mean, sd, p50, p25, p75, min, max, n)
	}
	
	sum num_alt
	scalar num_alt_fired = r(sum)
	sum np
	scalar num_alt_responded = r(sum)
	dis "number to fire (before being responded)", num_alt_fired/num_alt_responded
	
	collapse (count) num_enc = encounter_num ///
	(mean) time4alert_pp_mean = time4alert_pp time4alert_np_mean = time4alert_np ///
	time4alert_dis_mean = time4alert_dis time4alert_ref_mean = time4alert_ref ///
	time4alert_mean = time4alert ///
	num_alt_mean = num_alt num_rs_mean = np ///
	(sd) time4alert_pp_sd = time4alert_pp time4alert_sd = time4alert ///
	num_alt_sd = num_alt num_pp_sd = np ///
	(sum) time4alert_pp = time4alert_pp time4alert_np = time4alert_np time4alert = time4alert ///
	time4alert_dis = time4alert_dis ///
	num_alt = num_alt num_pp = pp num_np = np num_dis_ref = dis_or_ref ///
	num_enc_pp = enc_pp num_enc_np = enc_np num_enc_dis = enc_dis ///
	, by(clinic)
	
	gen num2fire = num_alt / num_dis_ref
	gen pp_time_rate = time4alert_pp / time4alert
	
	list
	keep clinic num2fire pp_time_rate num_alt_mean num_alt_sd ///
	time4alert_pp time4alert_np time4alert_dis time4alert num_enc_pp num_enc_np num_enc_dis num_enc
	order clinic num2fire pp_time_rate num_alt_mean num_alt_sd ///
	time4alert_pp time4alert_np time4alert_dis time4alert num_enc_pp num_enc_np num_enc_dis num_enc
	
	export excel using "TT_results.xlsx", sheet("TT_burden_`enc_type'") sheetreplace firstrow(variables)	
end



/*
foreach enc_type in allEnc { 
	log using "logfiles/analyze_TS_completion_`enc_type'.log", replace
	analyze_TS_completion `enc_type'
	log close
}
*/

/*
log using "logfiles/analyze_TT_completion_allEnc.log", replace
analyze_TT_completion
log close
*/


//*
foreach enc_type in allEnc { 
	log using "logfiles/analyze_TS_burden_`enc_type'.log", replace
	analyze_TS_burden `enc_type'
	log close
}

log using "logfiles/analyze_TT_burden_allEnc.log", replace
analyze_TT_burden
log close
*/
