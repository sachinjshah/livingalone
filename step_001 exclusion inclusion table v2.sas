
/* ---------------------------------------------------------

create tables for cohort flow

------------------------------------------------------------*/

/* Call in data */
OPTIONS nofmterr;
libname in "C:\Users\sachi\Box Sync\HRS\processed data";
%let cont_vars = 	AGE ;
%let cat_vars =  	RAFEMALE cog /*RAHISPAN*/  RARACEM HSGRAD INC_ABV_MED
							poorhealth PAIN BADVISION  BADHEARING
  							DIABETES HIGHBP  CANCER  LUNGD  HEARTD  STROKE
  							ARTHRITIS  DEPRESSION  
 							ADLdiff /*PROXY*/;
PROC FORMAT;
VALUE  eli
	1="Met inclusion criteria"
	0.1="Resides in NH (did not meet inclusion criteria)"
	0.2="< 65yrs (did not meet inclusion criteria)"
	0.3="Get's help with and ADL (did not meet inclusion criteria)"
	0.4="Missing ADL data (excluded)"
	0.5="Get's help with an IADL or other activity (did not meet inclusion criteria)"
	0.6="Did not answer support question (excluded)"
	0.7="Proxy interview (excluded)"
	0.8="Missing covariate data (excluded)";

run;

*exclusion;
data a; set in.analytic_long_2020_05_08;
where LA = "Yes";
if G097 in (5, 8) then nosupport = 1;
if G097 = 1 then nosupport = 0;
shock = min(shock, 1);

if eli = 1 and proxy = 1 then eli = 0.7;

array vv RAFEMALE RARACEM RAHISPAN married
		HSGRAD INC_ABV_MED poorhealth
		PAIN BADVISION  BADHEARING
		HIGHBP DIABETES CANCER  LUNGD  
		HEARTD  STROKE ARTHRITIS  DEPRESSION  
		SMOKING 
		BEDDIFF BATHDIFF DRESSDIFF EATDIFF
		TOILTDIFF WALKDIFF
		COG PROXY;
do over vv;
	if vv = . and eli = 1 then eli = 0.8;
end;
format eli eli.;
run;

proc sort data = a;
by wave;
run;

ods output OneWayFreqs = d;
proc freq data = a;
tables eli;
by wave;
where (7<wave<13) and eli in (1, 0.4, 0.6, 0.7, 0.8);
run;

data d; set d;
drop CumPercent Percent table F_ELI CumFrequency;
run;

proc sort data = d;
by eli;
run;

proc transpose data=d out=d_wide prefix=Wave;
    by eli ;
    id wave;
    var Frequency;
run;

data d_wide;
set d_wide;
drop _NAME_;
array ww wave:;
do over ww;
	if ww = . then ww = 0;
end;
run;

PROC EXPORT DATA= WORK.d_wide 
            OUTFILE= "C:\Users\sachi\Box Sync\HRS\results\cohort exclusion 2020-05-20.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

ods output OneWayFreqs = e;
proc freq data = a;
tables out:  /missing;
by wave;
where (7<wave<13) and eli = 1;
run;


data e; set e;
length Outcome $14.;
if OUT_ADLhelp ne .C and OUT_deceased ne .C and OUT_nursinghome ne .C then delete;
if OUT_ADLhelp = .C then outcome = "ADL dependency";
if OUT_deceased = .C then outcome = "Deceased";
if OUT_nursinghome = .C then outcome = "Nursing home";
drop CumPercent Percent table F_ELI CumFrequency F_: out_: ;
run;

proc sort  data = e; 
by outcome;
run;

proc transpose data=e out=e_wide prefix=Wave;
    by outcome ;
    id wave;
    var Frequency;
run;

PROC EXPORT DATA= WORK.e_wide 
            OUTFILE= "C:\Users\sachi\Box Sync\HRS\results\cohort exclusion outcome 2020-05-13.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


*save exclusion as a data set;
data a1; set a;
where eli = 1 and LA = "Yes" and (7<wave<13);
run;

data in.analytic_long_2020_05_08_no_miss;
set work.a1;
run;
