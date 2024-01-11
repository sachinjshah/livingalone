/* ---------------------------------------------------------

Program created by Sachin J Shah

v2 program goals: 
* run adjusted model 
* 3 outcomes (ADL, NH, died) 
* just looking at living alone
* calculate marginal effects

------------------------------------------------------------*/

/* Call in data */
OPTIONS nofmterr;
libname in "M:\HRS\processed data";
%let cont_vars = 	AGE ;
%let cat_vars =  	RAFEMALE cog /*RAHISPAN*/  RARACEM HSGRAD INC_ABV_MED
							poorhealth PAIN BADVISION  BADHEARING
  							DIABETES HIGHBP  CANCER  LUNGD  HEARTD  STROKE
  							ARTHRITIS  DEPRESSION  
 							ADLdiff /*PROXY*/;

							proc freq data = a; tables RAFEMALE;run;

data a; set in.long_20200808_NHsensi_nomiss;
run;
%INCLUDE "M:\!Macros\MARGINS_ver108.sas";


%MACRO adj_outs(out, out_name);
%Margins(	data      	= a,
            response  	= &out,
			class		= nosupport hhidpn &cat_vars,
            model     	= nosupport &cont_vars &cat_vars,
            dist      	= binomial,
			link 	  	= logit,
			margins    	= nosupport,
			geesubject	= hhidpn,
			geecorr		= cs,
            options   	= cl desc diff nomodel)

data marg_&out_name; 
set _margins;
length outcome $4.;
length model $11.;
adjusted = "Y";
outcome = "&out_name";
run;

data diffs_&out_name; 
set _diffs;
length outcome $4.;
length model $11.;
adjusted = "Y";
outcome = "&out_name";
run;

data marg_&out_name;
set marg_&out_name;
if outcome = "" then delete;
run;

data diffs_&out_name;
set diffs_&out_name;
if outcome = "" then delete;
run;

proc datasets nolist;
delete _:;
run;
quit;
%mend;

%adj_outs(OUT_nursinghome, NH);
%adj_outs(OUT_deceased, Died);
%adj_outs(OUT_ADLhelp, ADL);


data _marg; set marg_:;
if nosupport = 1 then support = "No support";
if nosupport = 0 then support = "Support";
array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.001);
end;
drop i;
output = cat(estimate*100, "% (", lower*100, "% to ",upper*100, "%)");
if outcome = "ADL" then n = 1;
else if outcome = "NH" then n = 2;
else n=3;
keep support outcome n adjusted output;
run;

proc sort; by n support; run;
proc print data = _marg; run;

data _diffs; set diffs_:;
array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.001);
end;
drop i;
output = cat(diff*-100, "% (", upper*-100, "% to ",lower*-100, "%)");
if outcome = "ADL" then n = 1;
else if outcome = "NH" then n = 2;
else n=3;
keep outcome n adjusted output;
run;

proc sort; by n; run;
proc print; run;



proc datasets nolist;
delete Table_: marg_: diffs_: clean: ;
run;
quit;



/* save coefficienets*/

%macro coeff(out, name, shock);
proc genmod data = a descending  ;
	class hhidpn &cat_vars /Ref = first;
	model &out = nosupport &shock &cat_vars &cont_vars/ dist=bin link=logit; 
	repeated subject=hhidpn / type=cs ;
	ods output GEEEmpPEst=b1;
run;

data b2;
length Parameter $40.;
length Level $30.;
set b1;
if z = . then delete;
Parameter = parm;
if Parm = "nosupport" then Parameter = "No support";
else if Parm = "SHOCK" then Parameter = "Health shock";
else if Parm = "nosupport*SHOCK" then Parameter = "No support by Health shock interaction";
else if Parm = "RAFEMALE" then Parameter = "Female";
else if Parm = "COG" then do;
	Parameter = "Cognition";
	if level1 = 2 then level = "Impairment not dementia";
	if level1 = 3 then level = "Dementia";
	end;
else if Parm = "RARACEM" then do;
	Parameter = "Race";
	if level1 = 2 then level = "Black";
	if level1 = 3 then level = "Other";
	end;
else if Parm = "HSGRAD" then Parameter = "High school graduate";
else if Parm = "INC_ABV_MED" then Parameter = "Net worth more than median";
else if Parm = "poorhealth" then Parameter = "Fair or poor self-reported health";
else if Parm = "PAIN" then Parameter = "Significant pain";
else if Parm = "BADVISION" then Parameter = "Visual impairment";
else if Parm = "BADHEARING" then Parameter = "Hearing impairment";
else if Parm = "DIABETES" then Parameter = "Diabetes";
else if Parm = "HIGHBP" then Parameter = "Hypertension";
else if Parm = "CANCER" then Parameter = "Cancer";
else if Parm = "LUNGD" then Parameter = "Lung disease";
else if Parm = "HEARTD" then Parameter = "Heart disease";
else if Parm = "STROKE" then Parameter = "Stroke";
else if Parm = "ARTHRITIS" then Parameter = "Arthritis";
else if Parm = "DEPRESSION" then Parameter = "Depression";
else if Parm = "ADLDIFF" then Parameter = "Difficulty with any ADL";
else if Parm = "AGE" then Parameter = "Age";
else if Parm = "" then Parameter = "";
drop parm level1 Stderr Z;
label Estimate = "Parameter estimate"
		Stderr = "Parameter standard error";
run;

proc print data = b2 label noobs;
title "&name";
run;

%mend;

%coeff(OUT_ADLhelp, ADL dependence, );
%coeff(out_nursinghome, Nursing home stay, );
%coeff(OUT_deceased, Mortality, );
%coeff(OUT_ADLhelp, ADL dependence w/ health shock interaction analysis, shock nosupport*shock);
%coeff(out_nursinghome, Nursing home stay w/ health shock interaction analysis, shock nosupport*shock);
%coeff(OUT_deceased, Mortality w/ health shock interaction analysis, shock nosupport*shock);
