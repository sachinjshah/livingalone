
/* ---------------------------------------------------------

Program created by Sachin J Shah

program goals: 
* run unadjsuted model 
* 3 outcomes (ADL, NH, died) 
* using 2 different approaches -> log poisson model and logit
* just looking at living alone
* calculate marginal effects

v3 is an updated analysi after the coding error in the development file was fixed 8-8-21

------------------------------------------------------------*/

/* Call in data */
OPTIONS nofmterr;
libname in "C:\Users\sachi\Box Sync\HRS\processed data";
%let cont_vars = 	AGE ;
%let cat_vars =  	RAFEMALE cog /*RAHISPAN*/  RARACEM HSGRAD INC_ABV_MED
							poorhealth PAIN BADVISION  BADHEARING
  							DIABETES HIGHBP  CANCER  LUNGD  HEARTD  STROKE
  							ARTHRITIS  DEPRESSION  
							ADLdiff
							/*PROXY*/;

data a; set in.analytic_long_20200808_nomiss;
/*data a; set in.analytic_long_2020_05_08_no_miss;*/
run;

%INCLUDE "C:\Users\sachi\Box Sync\data\!Macros\MARGINS.sas";

%MACRO adj_outs_shock(out, out_name);
ods text = "shock stratified outcomes for &out_name";

data marg_&out_name;
run;

data diffs_&out_name;
run;

%Margins(	data      	= a,
            response  	= &out,
			class		= nosupport shock hhidpn &cat_vars,
            model     	= nosupport shock nosupport*shock &cont_vars &cat_vars,
			dist      	= binomial,
			link 	  	= logit,
			margins    	= nosupport shock,
			geesubject	= hhidpn,
			geecorr		= cs,
            options   	= cl desc diff nomodel)

data _margins; 
set _margins;
length outcome $4.;
length model $11.;
adjusted = "Y";
model = "Logistic";
outcome = "&out_name";
run;

data _diffs; 
set _diffs;
length outcome $4.;
length model $11.;
adjusted = "Y";
model = "Logistic";
outcome = "&out_name";
run;

data marg_&out_name;
set _margins marg_&out_name;
if outcome = "" then delete;
run;

data diffs_&out_name;
set _diffs diffs_&out_name;
if outcome = "" then delete;
run;

proc datasets nolist;
delete _:;
run;
quit;
%mend;

%adj_outs_shock(OUT_nursinghome, NH);
%adj_outs_shock(OUT_deceased, Died);
%adj_outs_shock(OUT_ADLhelp, ADL);


data marg;
set marg_:;
run;

data diffs;
set diffs_:;
run;

data marg;
set marg;
format mu ll ul 5.1;
length support $10.;
length health_shock $15.;
where model = "Logistic";

if nosupport = 1 then Support = "No support";
else Support = "Support";

if shock = 0 then health_shock = "No health shock";
else health_shock = "Health shock";

array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.001);
end;
drop i;
Est_CI = cat(_mu*100, "% (", lower*100, "% to ",upper*100, "%)");
keep _mu lower upper Support health_shock outcome Est_CI ;
run;


data diffs;
set diffs;
where model = "Logistic" and Comp in ("1 - 3", "2 - 4");
if comp = "1 - 3" then comparison = "No health shock";
if comp = "2 - 4" then comparison = "Health shock";
array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.001);
end;
drop i;
Diff_CI = cat(diff*-100, "% (", upper*-100, "% to ",lower*-100, "%)");
Support = "Support";
p_value = Pr;
keep  comparison diff_CI support outcome p_value;
run;

proc print data = diffs; run;
proc print data = marg;run;

proc sql;
create table fin as 
select * from marg
left join diffs
on marg.outcome = diffs.outcome
and marg.health_shock = diffs.comparison 
and marg.support = diffs.support;

data fin;
length outcome $20.;
set fin;
if outcome = "ADL" then do; n = 1; outcome = "ADL dependency"; end;
else if outcome = "NH" then do; n = 2; outcome = "Nursing home stay"; end;
else do; n = 3; outcome = "Deceased"; end;
drop comparison;
run;

proc sort;
by n descending health_shock  support;
run;

data fin; set fin;
drop n;
label est_ci = "Predicted likelihood (95% CI)"
	  Diff_ci = "Average marginal effect (95% CI)"
	  p_value = "P value"
	health_shock = "Health shock"
	support = "Support"
	outcome = "Outcome";
run;

proc print data = fin label noobs;
run;

PROC EXPORT DATA= WORK.fin 
            OUTFILE= "C:\Users\sachi\Box Sync\HRS\results\Adjusted outcomes support x shock.csv" 
            DBMS=CSV label REPLACE;
     PUTNAMES=YES;
RUN;

proc datasets nolist;
delete Table_: marg: diffs: clean: ;
run;
quit;





