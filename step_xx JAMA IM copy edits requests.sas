/* ---------------------------------------------------------

Program created by Sachin J Shah
5/12/20

v2 program goals: 
* run unadjsuted model 
* 3 outcomes (ADL, NH, died) 
* just looking at living alone
* calculate marginal effects

------------------------------------------------------------*/

/* Call in data */
OPTIONS nofmterr;
libname in "N:\HRS\processed data";

data a; set in.long_20200808_NHsensi_nomiss; run;

proc contents; run;

proc print data = a (obs = 10);run;

proc sort data = a;
by hhidpn wave;
run;

data first; set a;
by hhidpn wave;
if first.hhidpn = 1;
run;

proc means data = first median p25 p75 mean n ;
var age rafemale;
run;

proc freq data = first;
tables rafemale * support;
run;

proc freq data = a;
table shock / missing;
run;

proc freq data = a;
tables hhidpn * shock / out = test noprint;
run;

proc print data = test (obs = 10); run;

proc freq data = test nlevels;
tables hhidpn / noprint;
run;

proc freq data =test;
tables shock;
run;

