
/* ---------------------------------------------------------

Program created by Sachin J Shah

* goal to perform RERI analysis
* VanderWeele Epidemiol. Methods 2014; 3(1): 33–72
* Correia Am J Epidemiol. 2018;187(11):2470–2480

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


data a; set in.analytic_long_2020_05_08_no_miss;
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
UNADJUSTED RERI model 
----------------------------------*/
** Approach 1: Use the VanderWeele OR code NOT accouting for clustering;
proc logistic descending data = a outest = b covout;
	model &out = nosupport shock nosupport*shock;
	*ods output ConvergenceStatus = conv;
	ods output ParameterEstimates = b2;
run;

data rerioutput;
set b;
array mm{*}_numeric_;
b0 = lag4(mm[1]);
b1 = lag4(mm[2]);
b2 = lag4(mm[3]);
b3 = lag4(mm[4]);
v11 = lag2(mm[2]);
v12 = lag(mm[2]);
v13 = mm[2];
v22 = lag(mm[3]);
v23 = mm[3];
v33 =  mm[4];
k1 = exp(b1 + b2+ b3)-exp(b1);
k2 = exp(b1+ b2+ b3)-exp(b2);
k3 =  exp(b1+ b2+ b3);
vreri =  v11*k1*k1+ v22*k2*k2+ v33*k3*k3+ 2*v12*k1*k2+ 2*v13*k1*k3
+ 2*v23*k2*k3;
reri = exp(b1+ b2+ b3)-exp(b1)-exp(b2)+ 1;
se_reri = sqrt(vreri);
ci95_l = reri-1.96*se_reri;
ci95_u = reri + 1.96*se_reri;
keep reri se_reri ci95_l ci95_u;
if _n_ = 5;
run;

/* proc print data = rerioutput;
var reri se_reri ci95_l ci95_u;
run; */

*reri;
data rerioutput; set rerioutput;
Adjusted = "N";
Clustering = "N";
Model = "Logistic";
Outcome = "&out_name";
keep reri se_reri ci95_l ci95_u Model Clustering Adjusted Outcome;
run;

data table_&out_name.;
set table_&out_name. rerioutput;
run;

*risk on multiplicative scale;
data b2; set b2;
if variable = "nosupport*SHOCK";
Outcome = "&out_name";
Model = "logit";
Adjusted = "N";
Clustering = "N";
run;

data table_&out_name._multi; 
set table_&out_name._multi b2; 
run;

proc datasets nolist;
delete b b2 rerioutput conv;
run;
quit;

** Approach 2: GENMOD LOGISTIC;
proc genmod data = a descending  ;
	class hhidpn;
	model &out = nosupport shock nosupport*shock / dist=bin link=logit; 
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
Adjusted = "N";
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
Outcome = "&out_name";
Adjusted = "N";
Clustering = "Y";
Model = "GEE binary logit";
drop df WaldChiSq _ESTTYPE_ Level1 LowerCL UpperCL Z;
run;

data table_&out_name._multi;
set table_&out_name._multi b1;
run;

proc datasets nolist;
delete b1 c1 d1 rerioutput1 conv;
run;
quit;

** Approach 3: GENMOD MODIFIED POISSON ;

proc genmod data = a descending  ;
	class hhidpn;
	model &out = nosupport shock nosupport*shock / dist=Poisson link=log; 
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
vreri =  v11*k1*k1+ v22*k2*k2+ v33*k3*k3+ 2*v12*k1*k2+ 2*v13*k1*k3
+ 2*v23*k2*k3;
reri = exp(b1+ b2+ b3)-exp(b1)-exp(b2)+ 1;
se_reri = sqrt(vreri);
ci95_l = reri-1.96*se_reri;
ci95_u = reri + 1.96*se_reri;
*keep reri se_reri ci95_l ci95_u;
if _n_ = 4;
run;

data rerioutput1; set rerioutput1;
Adjusted = "N";
Clustering = "Y";
Model = "GEE log poisson";
Outcome = "&out_name";
keep reri se_reri ci95_l ci95_u Model Clustering Adjusted Outcome;
run;

data table_&out_name.;
set table_&out_name. rerioutput1;
run;

** excess risk on multiplicive scale;
data b1; set b1;
if Parm = "nosupport*SHOCK";
Adjusted = "N";
Clustering = "Y";
Model = "GEE log poisson";
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


** Approach 4: GENMOD LOG BINOMIAL;

proc genmod data = a descending  ;
	class hhidpn;
	model &out = nosupport shock nosupport*shock / dist=bin link=log; 
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
vreri =  v11*k1*k1+ v22*k2*k2+ v33*k3*k3+ 2*v12*k1*k2+ 2*v13*k1*k3
+ 2*v23*k2*k3;
reri = exp(b1+ b2+ b3)-exp(b1)-exp(b2)+ 1;
se_reri = sqrt(vreri);
ci95_l = reri-1.96*se_reri;
ci95_u = reri + 1.96*se_reri;
*keep reri se_reri ci95_l ci95_u;
if _n_ = 4;
run;

data rerioutput1; set rerioutput1;
Adjusted = "N";
Clustering = "Y";
Model = "GEE log binomial";
Outcome = "&out_name";
keep reri se_reri ci95_l ci95_u Model Clustering Adjusted Outcome;
run;

data table_&out_name.;
set table_&out_name. rerioutput1;
run;

** excess risk on multiplicive scale;
data b1; set b1;
if Parm = "nosupport*SHOCK";
Adjusted = "N";
Clustering = "Y";
Model = "GEE log binomial";
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

/*--------------------------------
ADJUSTED RERI model 
----------------------------------*/
** Approach 5: Use the VanderWeele OR code NOT accouting for clustering;
proc logistic descending data = a outest = b covout;
	class &cat_vars;
	model &out = nosupport shock nosupport*shock &cat_vars &cont_vars;
	ods output ParameterEstimates = b2;
run;

data rerioutput;
set b;
array mm{*}_numeric_;
b0 = lag4(mm[1]);
b1 = lag4(mm[2]);
b2 = lag4(mm[3]);
b3 = lag4(mm[4]);
v11 = lag2(mm[2]);
v12 = lag(mm[2]);
v13 = mm[2];
v22 = lag(mm[3]);
v23 = mm[3];
v33 =  mm[4];
k1 = exp(b1 + b2+ b3)-exp(b1);
k2 = exp(b1+ b2+ b3)-exp(b2);
k3 =  exp(b1+ b2+ b3);
vreri =  v11*k1*k1+ v22*k2*k2+ v33*k3*k3+ 2*v12*k1*k2+ 2*v13*k1*k3
+ 2*v23*k2*k3;
reri = exp(b1+ b2+ b3)-exp(b1)-exp(b2)+ 1;
se_reri = sqrt(vreri);
ci95_l = reri-1.96*se_reri;
ci95_u = reri + 1.96*se_reri;
keep reri se_reri ci95_l ci95_u;
if _n_ = 5;
run;

data rerioutput; set rerioutput;
Adjusted = "Y";
Clustering = "N";
Model = "Logistic";
Outcome = "&out_name";
keep reri se_reri ci95_l ci95_u Model Clustering Adjusted Outcome;
run;

data table_&out_name.;
set table_&out_name. rerioutput;
run;

*risk on multiplicative scale;
data b2; set b2;
if variable = "nosupport*SHOCK";
Adjusted = "Y";
Clustering = "N";
Model = "Logistic";
Outcome = "&out_name";
run;

data table_&out_name._multi; set table_&out_name._multi b2; run;

proc datasets nolist;
delete b b2 rerioutput conv;
run;
quit;

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

** Approach 7: GENMOD MODIFIED POISSON ;

proc genmod data = a descending  ;
	class hhidpn &cat_vars;
	model &out = nosupport shock nosupport*shock &cat_vars &cont_vars / dist=Poisson link=log; 
	repeated subject=hhidpn / type=cs ecovb;
	ods output GEEEmpPEst=b1;
	ods output GEERCov = c1;
	ods output ConvergenceStatus = conv;
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
vreri =  v11*k1*k1+ v22*k2*k2+ v33*k3*k3+ 2*v12*k1*k2+ 2*v13*k1*k3
+ 2*v23*k2*k3;
reri = exp(b1+ b2+ b3)-exp(b1)-exp(b2)+ 1;
se_reri = sqrt(vreri);
ci95_l = reri-1.96*se_reri;
ci95_u = reri + 1.96*se_reri;
*keep reri se_reri ci95_l ci95_u;
if _n_ = 4;
run;

data rerioutput1; set rerioutput1;
Adjusted = "Y";
Clustering = "Y";
Model = "GEE log poisson";
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
Model = "GEE log poisson";
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


** Approach 8: GENMOD LOG BINOMIAL -> Did not converge;

proc genmod data = a descending  ;
	class hhidpn &cat_vars;
	model &out = nosupport shock nosupport*shock &cat_vars &cont_vars/ dist=bin link=log; 
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
vreri =  v11*k1*k1+ v22*k2*k2+ v33*k3*k3+ 2*v12*k1*k2+ 2*v13*k1*k3
+ 2*v23*k2*k3;
reri = exp(b1+ b2+ b3)-exp(b1)-exp(b2)+ 1;
se_reri = sqrt(vreri);
ci95_l = reri-1.96*se_reri;
ci95_u = reri + 1.96*se_reri;
*keep reri se_reri ci95_l ci95_u;
if _n_ = 4;
run;

data rerioutput1; set rerioutput1;
Adjusted = "Y";
Clustering = "Y";
Model = "GEE log binary";
Outcome = "&out_name";
keep reri se_reri ci95_l ci95_u Model Clustering Adjusted Outcome;
run;

data table_&out_name.;
set table_&out_name. rerioutput1;
run;

data b1; set b1;
if Parm = "nosupport*SHOCK";
Adjusted = "Y";
Clustering = "Y";
Model = "GEE log binary";
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
%RERI(OUT_deceased, Died);
%RERI(OUT_ADLhelp, ADL);


data table_all;
set table_adl table_died table_nh;
run;

PROC EXPORT DATA= WORK.table_all
            OUTFILE= "C:\Users\sachi\Box Sync\HRS\results\RERI_add_all.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


data table_all_multi;
set table_adl_multi table_died_multi table_nh_multi;
run;

PROC EXPORT DATA= WORK.table_all_multi 
            OUTFILE= "C:\Users\sachi\Box Sync\HRS\results\RERI_multi_all.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;








