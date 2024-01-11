
/* ---------------------------------------------------------

program by sachin shah
goal to perform NH sensitivity analysis 
looking at 30, 60 and 100 day time intervals. 

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

data a; 
set in.long_20200808_NHsensi_nomiss;
/*set in.analytic_long_20200808_nomiss;*/
run;

*%INCLUDE "C:\Users\sachi\Box Sync\data\!Macros\MARGINS_ver108.sas";
%INCLUDE "M:\!Macros\MARGINS_ver108.sas";

proc freq data = a;
tables OUT:;
run;

%MACRO adj_outs(out, out_name, label);
data marg_&out_name;
run;

data diffs_&out_name;
run;

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

data _margins; 
set _margins;
length outcome $36.;
length model $11.;
adjusted = "Y";
model = "Logistic";
outcome = "&label";
run;

data _diffs; 
set _diffs;
length outcome $36.;
length model $11.;
adjusted = "Y";
model = "Logistic";
outcome = "&label";
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

%adj_outs(OUT_nursinghome_any, anyNH, Any nursing home days);
%adj_outs(OUT_nursinghome, NH, Nursing home 30 days (original));
%adj_outs(OUT_nursinghome60, NH60, Nusring home 60 days (sensitivity));
%adj_outs(OUT_nursinghome100, NH100, Nursing home 100 days (sensitivity));


data marg;
length support $10.;

set marg_nh marg_nh60 marg_nh100;
if nosupport = 1 then support = "No support";
if nosupport = 0 then support = "Support";

array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.001);
end;
drop i;
a = cat(estimate*100, "% (", lower*100, "% to ",upper*100, "%)");
keep Support outcome a n;
n+1;
run;


data diffs; 
set diffs_nh diffs_nh60 diffs_nh100;

array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.001);
end;
drop i;
Difference = cat(diff*-100, "% (", upper*-100, "% to ",lower*-100, "%)");
Support = "No support";
keep Support outcome Difference;

run;

proc print data = diffs;
run;

proc sql;
create table fin
as select* from marg
left join diffs
on marg.outcome = diffs.outcome and
marg.support = diffs.support;
quit;

proc sort data = fin;
by n;
run;

data fin; set fin;
rename support = Support
		Outcome = Outcome
		a = Estimate;
drop n;
run;
proc print data = fin;run;


proc datasets nolist lib = work;
delete marg: diff:;
run;
quit;

/*INTERACTION analysis*/

%MACRO adj_outs_shock(out, out_name, label);

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
/*			at		= cog */

data _margins; 
set _margins;
length outcome $35.;
length model $11.;
outcome = "&label";
run;

data _diffs; 
set _diffs;
length outcome $35.;
length model $11.;
outcome = "&label";
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

%adj_outs_shock(OUT_nursinghome_any, NHany, Any nursing home days);
%adj_outs_shock(OUT_nursinghome, NH, Nursing home 30 days (original));
%adj_outs_shock(OUT_nursinghome60, NH60, Nusring home 60 days (sensitivity));
%adj_outs_shock(OUT_nursinghome100, NH100, Nursing home 100 days (sensitivity));


data marg;
set marg_nhany marg_nh marg_nh60 marg_nh100;

if nosupport = 1 then Support = "No support";
else Support = "Support";

if shock = 0 then health_shock = "No health shock";
else health_shock = "Health shock";

array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.001);
end;
drop i;
Est = cat(_mu*100, "% (", lower*100, "% to ",upper*100, "%)");
keep Support health_shock outcome Est n;
if outcome = "Any nursing home days" then n = 1;
else if outcome = "Nursing home 30 days (original)" then n = 2;
else if outcome = "Nusring home 60 days (sensitivity)" then n = 3;
else n = 4;
run;

proc sort data = marg; by n health_shock;
run;

proc print data = marg;run;

data diffs;
set diffs_:;
run;

data diffs; 
length comparison $15.;
set diffs;
if comp = "1 - 3" then comparison = "No health shock";
if comp = "2 - 4" then comparison = "Health shock";
array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.001);
end;
drop i;
Difference = cat(diff*-100, "% (", upper*-100, "% to ",lower*-100, "%)");

keep outcome comparison Difference Support;
if comparison = "" then delete;
Support = "No support";
run;

proc sql;
create table fin_shock
as select * from marg
left join diffs
on marg.outcome = diffs.outcome 
and marg.support = diffs.support 
and marg.health_shock = diffs.comparison;
quit;

proc sort data = fin_shock;
by n  descending  health_shock  support;
run;

data fin_shock ; set fin_shock;
drop n comparison;
label health_shock = "Health Shock"
		outcome = "Outcome"
		est = "Estimate";
run;

proc print data = fin label; run;
proc print data = fin_shock label;run;


proc datasets nolist lib = work;
delete marg diff:;
run;
quit;



%macro RERI(out, out_name, data);

data table_&out_name.;
length Outcome $4;
length Model $16.;
length Adjusted $1.;
length Clustering $1.;
format RERI BEST12.;
format SE_RERI BEST12.;
format ci95_l BEST12.;
format ci95_u BEST12.;
run;

data table_&out_name._multi;
length Outcome $4;
length Model $16.;
length Adjusted $1.;
length Clustering $1.;
length variable $15.;
run;

/*--------------------------------
ADJUSTED RERI model 
----------------------------------*/

** Approach 6: GENMOD LOGISTIC;
proc genmod data = &data. descending  ; **data = b for ADL diff subset;
	class hhidpn &cat_vars;
	model &out = nosupport shock nosupport*shock &cat_vars &cont_vars/ dist=bin link=logit; 
	repeated subject=hhidpn / type=cs ecovb;
	ods output GEEEmpPEst=b1;
	ods output GEERCov = c1;
run;

data b1; set b1; id = _N_ ; run;
data c1; set c1; id = _N_; run;

proc sql;
create table d1 as 
select * from b1 
left join c1 
on b1.id = c1.id;
quit;

data rerioutput1;
set d1;
array mm{*}_numeric_;
b0 = lag3(mm[1]);
b1 = lag2(mm[1]);
b2 = lag1(mm[1]);
b3 = (mm[1]);
v11 = lag2(mm[9]);
v12 = lag(mm[9]);
v13 = (mm[9]);
v22 = lag(mm[10]);
v23 = (mm[10]);
v33 = (mm[11]);
k1 = exp(b1 + b2+ b3)-exp(b1);
k2 = exp(b1+ b2+ b3)-exp(b2);
k3 = exp(b1+ b2+ b3);
vreri =  v11*k1*k1 + v22*k2*k2 + v33*k3*k3 + 2*v12*k1*k2 + 2*v13*k1*k3 + 2*v23*k2*k3;
reri = exp(b1+ b2+ b3) - exp(b1) - exp(b2) + 1;
se_reri = sqrt(vreri);
ci95_l = reri - 1.96*se_reri;
ci95_u = reri + 1.96*se_reri;
*keep reri se_reri ci95_l ci95_u;
if _n_ = 4;
run;

data rerioutput1; set rerioutput1;
Adjusted = "Y";
Clustering = "Y";
Model = "GEE binary logit";
Outcome = "&out_name";
keep reri se_reri ci95_l ci95_u Model Clustering Adjusted Outcome;
run;

data table_&out_name.;
set table_&out_name. rerioutput1;
run;

** excess risk on multiplicive scale;
data b1; set b1;
if Parm = "nosupport*SHOCK";
Adjusted = "Y";
Clustering = "Y";
Model = "GEE binary logit";
Outcome = "&out_name";
drop df WaldChiSq _ESTTYPE_ Level1 LowerCL UpperCL Z;
run;

data table_&out_name._multi; 
set table_&out_name._multi b1;
run;

proc datasets nolist;
delete b1 c1 d1 rerioutput1 conv;
run;
quit;

**clean up data;
data table_&out_name.;
set table_&out_name.;
if outcome = "" then delete;
array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.1);
end;
drop i SE_RERI;
run;

data table_&out_name._multi; 
set table_&out_name._multi; 
drop variable parm df id ClassVal0 _ESTTYPE_ WaldChiSq; 
if outcome = "" then delete ;
OR = exp(estimate);
LL = exp(estimate -1.96*Stderr);
UL = exp(estimate +1.96*Stderr);
array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.1);
end;
drop i Stderr;
run;

proc print data = table_&out_name. noobs label ;run;
proc print data = table_&out_name._multi noobs label; run;

%mend RERI;

%RERI(out_nursinghome_any, NH, a);
%RERI(out_nursinghome, NH, a);
%RERI(out_nursinghome60, NH, a);
%RERI(out_nursinghome100, NH, a);

