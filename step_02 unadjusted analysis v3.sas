/* ---------------------------------------------------------

Program created by Sachin J Shah

program goals: 
* run unadjsuted model 
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


data a; set in.long_20200808_NHsensi_nomiss;
run;
%INCLUDE "M:\!Macros\MARGINS_ver108.sas";

%MACRO unadj_outs(out, out_name);

data marg_&out_name;
run;

data diffs_&out_name;
run;

%Margins(	data      	= a,
            response  	= &out,
			class		= nosupport hhidpn,
            model     	= nosupport,
            dist      	= binomial,
			link 	  	= logit,
			margins    	= nosupport,
			geesubject	= hhidpn,
			geecorr		= cs,
            options   	= cl desc diff nomodel)

data _margins; 
set _margins;
length outcome $4.;
length model $11.;
adjusted = "N";
outcome = "&out_name";
run;

data _diffs; 
set _diffs;
length outcome $4.;
length model $11.;
adjusted = "N";
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

%unadj_outs(OUT_nursinghome, NH);
%unadj_outs(OUT_deceased, Died);
%unadj_outs(OUT_ADLhelp, ADL);


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



