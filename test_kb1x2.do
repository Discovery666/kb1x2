/*--------------------------------------------------------*/
// set up
clear all
macro drop _all
global jtmp "H:\kb1x2\data\jtmp"
cd "H:\kb1x2"
quietly do "code\i888_kb1x2.do"
/*--------------------------------------------------------*/
// generate a sample data set for our experiment
{
	 clear
	 set obs 10000
	 //
	 set seed 123456
	 gen ri = _n
	 gen x1 = round(rnormal()*10)
	 gen x2 = round(rnormal()*10) 
	 gen x3 = round(rnormal()*10) 
	 //
	 gen z1 = round(rnormal()*2) + x1 + 0.5*x2 + 0.5*x3
	 //
	 gen y = 1+ 1*z1 + x1 + 0.5*x2 + 0.5*x3 + 0.5*rnormal()
	 //
	 save "${jtmp}/DecomFigureA.dta", replace
	 }
/*--------------------------------------------------------*/
// define all the global macros to be used in the decomposition
{
	global SA x1 // context specific
	global SB x2 x3 // context specific
	//
	global tmacro depvar x1all x2all ClusterV dropped
	global tmacroi 1 2 3 4 999
	//
	global depvar y // you dependent variable
	global x1all z1 // your variable of interest
	global x2all0 x1 x2 x3 // the unrestricted set of variables
	global x2deltaVarlist0 = "SA SB" // the set of variables for decomposition, context specific
	global ClusterV ""	// cluster variable for robust r2
	}
/*--------------------------------------------------------*/
// call the decomposition program
k_b1x2 1
/*
------------------------------------------------------------------------------
           y |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
          SA |   .6537361   .0050633   129.11   0.000     .6438123    .6636599
          SB |   .3230493   .0047836    67.53   0.000     .3136737    .3324249
         _TC |   .9767854   .0024527   398.25   0.000     .9719782    .9815926
------------------------------------------------------------------------------
depvar:          y
x1all:           z1
*/
/*--------------------------------------------------------*/
// inspect the estimated coefficients of the restricted and unrestricted regression
use "${jtmp}/DecomFigureA.dta", clear
reg y z1 // 1.974672
reg y z1 x1 x2 x3 // .9978862
/*--------------------------------------------------------*/