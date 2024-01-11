
/* ---------------------------------------------------------

updated sensitivity analysis to look at robustness of NH definition
------------------------------------------------------------*/

/*
Call in data
*/

OPTIONS nofmterr;
libname in "C:\Users\sachi\Box Sync\HRS\processed data";
%let cont_vars = 	AGE ;
%let cat_vars =  	RAFEMALE cog /*RAHISPAN*/  RARACEM HSGRAD INC_ABV_MED
							poorhealth PAIN BADVISION  BADHEARING
  							DIABETES HIGHBP  CANCER  LUNGD  HEARTD  STROKE
  							ARTHRITIS  DEPRESSION  
 							ADLdiff /*PROXY*/;


data a; 
/*set in.analytic_long_NH_sensi_no_miss;*/
set in.long_20200808_NHsensi_nomiss;
run;


/*
Basic frequency tables
*/
** basic outcomes;
proc freq;
tables nosupport * OUT: / nocol nopercent;
run;

**interaction outcomes;
proc freq;
tables shock *nosupport *  out: / nocol nopercent nofreq missing;
run;



%macro RERI(out, out_name);

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
proc genmod data = a descending  ;
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
run;

data table_&out_name._multi; 
set table_&out_name._multi; 
drop variable parm df id ClassVal0 _ESTTYPE_ WaldChiSq; 
if outcome = "" then delete ;
run;

%mend RERI;

%RERI(OUT_nursinghome, NH);
%RERI(OUT_nursinghome60, NH60);
%RERI(OUT_nursinghome100, NH100);


data table_all;
set table_nh table_nh60 table_nh100;
if outcome = "NH" then Description = "Nursing home 30 days (original)    ";
else if outcome = "NH60" then Description = "Nusring home 60 days (sensitivity)";
else Description = "Nusring home 100 days (sensitivity)";

array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.1);
end;
drop i se_reri reri ci: adjusted clustering outcome model;
Est = cat(RERI, " (", ci95_l, " to ",ci95_u, ")");
label est = "RERI OR (95% CI )";

run;

proc print data = table_all label;
run;

data table_all_multi;
set table_nh_multi table_nh60_multi table_nh100_multi;
if outcome = "NH" then Description = "Nursing home 30 days (original)    ";
else if outcome = "NH60" then Description = "Nusring home 60 days (sensitivity)";
else Description = "Nusring home 100 days (sensitivity)";
OR = exp(estimate);
OR_LL = exp(estimate - 1.96 * stderr);
OR_UL = exp(estimate +1.96 * stderr);

array _nums {*} _numeric_;
do i = 1 to dim(_nums);
  _nums{i} = round(_nums{i},.01);
end;
*drop i se_reri reri ci: adjusted clustering outcome model;
Est = cat(OR, " (", OR_LL, " to ",OR_UL, ")");
label est = "Multiplicative interaction OR (95% CI)";
keep Description est;
run;

proc print data = table_all_multi label noobs;
run;

