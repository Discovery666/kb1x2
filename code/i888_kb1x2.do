/*--------------------------------------------------------*/
* (21) k_datetime
{
capture program drop k_datetime 
program define k_datetime 
	local t = "$S_DATE $S_TIME"
	local CurrentTime = clock("`t'","DMYhms")
	local StartTime = clock("$StartTime","DMYhms")
	local TimeLapsed = `CurrentTime' - `StartTime'
	local thh = mod(floor(hours(`TimeLapsed')),24)
	local tmm = mod(floor(minutes(`TimeLapsed')),60)
	local tss = mod(round(seconds(`TimeLapsed'),0.1),60)
	disp "`t' --- `thh':`tmm':`tss' passed so far"
end
}
/*--------------------------------------------------------*/
// SaveMatrixList
{
	capture program drop SaveMatrixList
	program define SaveMatrixList
	args MatrixList FileNameBase
	foreach ti of global `MatrixList' {
		SaveMatrix `ti' `FileNameBase'
		} // end foreach ti
	end
	}
/*--------------------------------------------------------*/
// SaveMatrix
{
	capture program drop SaveMatrix
	program define SaveMatrix
	args tMatrix FileNameBase
	preserve
	{
		clear
		svmat `tMatrix'
		save "${jtmp}/`FileNameBase'`tMatrix'.dta", replace
		} 
	end
	}
/*--------------------------------------------------------*/	
// LoadMatrixList
{
	capture program drop LoadMatrixList
	program define LoadMatrixList
	args MatrixList FileNameBase
	foreach ti of global `MatrixList' {
		LoadMatrix `ti' `FileNameBase'
		} // end foreach ti
	end
	}
/*--------------------------------------------------------*/
// LoadMatrix
{
	capture program drop LoadMatrix
	program define LoadMatrix
	args tMatrix FileNameBase
	preserve
	use "${jtmp}/`FileNameBase'`tMatrix'.dta", clear
	mkmat *, matrix(`tMatrix')
	end
	}
/*--------------------------------------------------------*/
// k_b1x2 
{
	capture program drop k_b1x2 
	program define k_b1x2 
	args tSeq
	// section 7.4
	{
		global x2deltaVarlist "${x2deltaVarlist0}"
		global x2all = "${x2all0}"
		}
	/*--------------------------------------------------------*/
	// section 7.5: generate variable list data set
	{
			// DecomFigureA.dta => tvarlistA.dta
			$quietlyFlag {
				use "${jtmp}/DecomFigureA.dta", clear
				keep if _n<10
				xpose, clear varname
				keep _varname
				rename _varname vname
				gen vi =_n
				save "${jtmp}/tvarlistA.dta", replace
				} 
			/*--------------------------------------------------------*/
			// tvarlistA.dta => tvarlistB.dta
			$quietlyFlag { // tvarlistB start
				use "${jtmp}/tvarlistA.dta", clear
				expand 1+(vi==1)
				sort vi
				quietly replace vi=0 if vi==1 & _n==1
				quietly replace vname = "_cons" if vi==0
				//
				quietly gen byte vtype = .
				quietly gen byte dropped = 0
				quietly gen x2groupi = .
				quietly gen x2groupnm = ""
				gen gj = .
				//
				local tmacrocount: word count $tmacro
				local tj = 0
				while `tj'<`=`tmacrocount'-1' {
					local tj = `tj'+1
					local tmname: word `tj' of $tmacro
					local tvtype: word `tj' of $tmacroi
					disp "`tmname':`tvtype'"
					if "`tmname'"~="dropped" {
						foreach tv of global `tmname' {
							quietly replace vtype = `tvtype' if vname==`"`tv'"'
							} // foreach tv end
						local tstr = usubstr("${`tmname'}",1,200)
						disp "`tmname': `tstr'..."
						}
					} // while tj end
				quietly replace vtype = 2 if vname=="_cons"
				quietly replace vtype= 999 if missing(vtype) 
				quietly replace dropped = 1 if vtype==999
				sort vtype vi
				/*--------------------------------------------------------*/
			$quietlyFlag { // x2 delta varlist setup start
					local groupi = 0
					foreach tv of global x2deltaVarlist {
						disp `" `tv' "', _continue
						$quietlyFlag {
							local groupi = `groupi'+1
							local jj = 0
							quietly foreach tw of global `tv' {
								local jj = `jj'+1
								replace x2groupnm = `"`tv'"' if vname==`"`tw'"' & vtype==3
								replace x2groupi = `groupi' if vname==`"`tw'"' & vtype==3
								replace gj = `jj' if vname==`"`tw'"' & vtype==3
								} // end foreach tw
							} // end quietly
						} // end foreach tv
					quietly replace x2groupnm = "Residual" if missing(gj) & vtype==3
					$quietlyFlag levelsof x2groupnm, local(tx2groupnm) clean
					global groupN: word count `tx2groupnm'
					quietly replace x2groupi = $groupN if x2groupnm == "Residual" & missing(gj) & vtype==3
					//
					sort vtype x2groupi	gj vi
					quietly by vtype x2groupi: replace gj = _n if vtype==3
					} // x2 delta varlist setup end
				/*--------------------------------------------------------*/
			 $quietlyFlag { // final organize start
					sort vtype x2groupi	gj vi
					quietly replace vi = _n
					//
					order vi vtype dropped x2groupi x2groupnm gj vname
					label define vtypelb 1 "dependent variable" 2 "X1" 3 "X2" 4 "cluster" 999 "not used", modify
					label value vtype vtypelb
					label define droppedlb 1 "dropped" 0 "use", modify
					label value dropped droppedlb
					format vi x2groupi gj %-4.0f
					format vtype %-20.0f
					format dropped %-12.0f
					format x2groupnm vname %-20s
					quietly compress
					} // final organize end
				save "${jtmp}/tvarlistA.dta", replace
				} // tvarlistB end
			/*--------------------------------------------------------*/
			// update the b1xa related variable lists
			$quietlyFlag {	// update variable lists start
				use "${jtmp}/tvarlistA.dta", clear
				local tj = 0
				while `tj'<5 {
					local tj = `tj'+1
					local tmname: word `tj' of $tmacro
					local tvtype: word `tj' of $tmacroi
					disp "`tmname':`tvtype'"
					//
					{
						use "${jtmp}/tvarlistA.dta", clear
						quietly keep if vtype==`tvtype' & vname~="_cons"
						keep vname vi
						sort vi
						local tempvar = ""
						local ti 0
						local timax = _N
						while `ti'< `timax' {
							local ti = `ti' + 1
							local tvname = vname[`ti']
							local tempvar = "`tempvar' `tvname' "			
							} // end while ti
						} // end this block
						capture global `tmname' = "`tempvar'"
					} // end while tj		
				} // update variable lists end
			capture macro list dropped
			/*--------------------------------------------------------*/
			// block 2: group list: groupnamesAeq groupnamesBeq
			{
				clear
				global groupnamesAeq = ""
				local tgi = 0
				while `tgi' < $groupN {
					local tgi = `tgi'+1
					global groupnamesA`tgi' = ""
					global groupnamesAeq`tgi' = ""
					//
					use "${jtmp}/tvarlistA.dta", clear
					keep vtype vname x2groupnm x2groupi gj 
					quietly keep if vtype==3 & x2groupi==`tgi'
					local tgmacro = "groupnamesA`tgi'"
					local tgEqmacro = "groupnamesAeq`tgi'"
					//
					sort x2groupi gj
					gen vnameeq = x2groupnm+":"+vname
					gen ri=_n
					local timax = _N
					local ti = 0
					while `ti'<`timax' {
						local ti = `ti'+1
						local tvname = vname[`ti']
						local tveqname = vnameeq[`ti']
						global `tgmacro' = "$`tgmacro' `tvname' "
						global `tgEqmacro' = "$`tgEqmacro' `tveqname' "
						} // end while ti
					global groupnamesAeq = "${groupnamesAeq} $`tgEqmacro'"
					clear
					} // end while tgi
				//
				clear
				global groupnamesBeq = ""
				local tgi = 0
				while `tgi'<$groupN {
					local tgi = `tgi'+1
					global groupnamesBeq`tgi' = ""
					//
					use "${jtmp}/tvarlistA.dta", clear
					keep vtype vname x2groupnm x2groupi gj
					quietly keep if vtype==3 & x2groupi==`tgi' & gj==1
					local tgEqmacro = "groupnamesBeq`tgi'"
					//
					sort x2groupi gj
					local tveqname = x2groupnm+":"
					foreach tvi in $x1all "_cons" {
						local tvj = "`tveqname'`tvi'"
						global `tgEqmacro' = "$`tgEqmacro' `tvj'"
						} // end foreach tvi
					global groupnamesBeq = "$groupnamesBeq $`tgEqmacro'"
					clear
					} // end while tgi
				} // end this block 2
			/*--------------------------------------------------------*/
			// block 3: key parameters
			{
				global k1: word count $x1all _con
				global k2: word count $x2all
				global k1S = 2
				global k1E = $k1+1
				global k2S = $k1+2
				global k2E = $k1+$k2+1		
				global k2Sb = $k1+1
				global k2Eb = $k1+$k2
				} // block 3 end
			/*--------------------------------------------------------*/
			// create the updated main data file DecomFigureB.dta based on updated variable lists
			$quietlyFlag {
				use "${jtmp}/DecomFigureA.dta", clear
				order $depvar $x1all $x2all $clusterV
				keep $depvar $x1all $x2all $clusterV
				save "${jtmp}/DecomFigureA.dta", replace	
				}
			/*--------------------------------------------------------*/
			disp "(1) initialization done"
			k_datetime
		} // 
	/*--------------------------------------------------------*/
	// section 7.6: calculate basic matrices
	{
			// block 4.1: calculate the matrics
			use "${jtmp}/DecomFigureA.dta", clear
			// calculate all the matrics
			{
				gen byte cons = 1
				matrix accum XYXY = $depvar $x1all cons $x2all, noconstant
				//
				global sst = XYXY[1,1]
				//
				matrix XX = XYXY[$k1S..$k2E,$k1S..$k2E]
				matrix XY = XYXY[$k1S..$k2E,1]
				matrix X1X1 = XYXY[$k1S..$k1E,$k1S..$k1E]
				matrix X1Y = XYXY[$k1S..$k1E,1]
				matrix X1X2 = XYXY[$k1S..$k1E,$k2S..$k2E]
				//
				matrix XXInv = syminv(XX)
				matrix VU = XXInv
				matrix colnames VU = $x1all cons $x2all
				matrix rownames VU = $x1all cons $x2all
				//
				matrix X1X1Inv = syminv(X1X1)
				matrix VR = X1X1Inv
				matrix colnames VR = $x1all cons 
				matrix rownames VR = $x1all cons 
				//
				matrix BetaU = (XXInv*XY)'
				matrix colnames BetaU = $x1all cons $x2all
				matrix rownames BetaU = $depvar
				//
				matrix BetaUX1 = BetaU[1,1..$k1]
				matrix BetaUX2 = BetaU[1,$k2Sb..$k2Eb]
				//
				matrix BetaR = (X1X1Inv*X1Y)'
				matrix colnames BetaR = $x1all cons
				matrix rownames BetaR = $depvar
				//
				matrix Gamma = (X1X1Inv*X1X2)
				matrix colnames Gamma = $x2all
				matrix rownames Gamma = $x1all cons
				//
				matrix Delta = BetaUX2*Gamma'
				matrix colnames Gamma = $x1all cons
				matrix rownames Gamma = $depvar
				} // end of block 4.1
			} // end step 4
	/*--------------------------------------------------------*/
	// section 7.7: calculate variance-covariance matrices
	{
		// block 5.2
		$quietlyFlag {
			use "${jtmp}/DecomFigureA.dta", clear
			gen byte cons =1 
			matrix score double yhat = BetaU
			gen double ehat = $depvar - yhat
			matrix score double yhatR = BetaR
			gen double ehatR = $depvar - yhatR
			save "${jtmp}/DecomFigureA.dta", replace
			} // end block 5.2
		// block 5.3 
		if "$clusterV"=="" {
			use "${jtmp}/DecomFigureA.dta", clear
			egen t = sum(ehat^2)
			global sse = t[1]
			global s2U = $sse /(_N-$k1-$k2)
			matrix VU = $s2U * XXInv
			drop t
			//
			egen t = sum(ehatR^2)
			global sseR = t[1]
			global s2R = $sseR /(_N-$k1)
			matrix VR = $s2R * X1X1Inv
			drop t
			//
			}
		else {
			use "${jtmp}/DecomFigureA.dta", clear
			capture gen byte cons=1
			quietly _robust ehat, v(VU) minus(`=$k1+$k2') cluster($clusterV)
			global NCluster = r(N_clust)
			//
			quietly _robust ehatR, v(VR) minus($k1) cluster($clusterV)
			} // end block 5.3
		// block 5.4
		{
			matrix VU1 = VU[1..$k1, 1..$k1] 
			matrix VU2 = VU[$k2Sb..$k2Eb, $k2Sb..$k2Eb]
			} // end block 5.2
		//
		disp "(2) variance-covariance matrices produced"
		k_datetime
		}
	/*--------------------------------------------------------*/
	// section 7.8: calculate Hand Calculation of biases variance-covariance
	{
		//global matrixList XXBigInv XXBig MVlong MV BigMat  AuxRegV CBigAux __TCcov _tcregpart __TCb
	/*--------------------------------------------------------*/
		// block 6.1: produce the BigMat matrix 
		{
			/*--------------------------------------------------------*/
			//// BigMat
			$quietlyFlag {
				clear
				matrix colnames Gamma = $x2all
				matrix rownames Gamma = $x1all cons
				quietly svmat Gamma, names(col)
				quietly xpose, varname clear
				rename _varname vname
				quietly merge 1:1 vname using "${jtmp}/tvarlistA.dta", keep(3) nogen keepusing(gj x2groupi )
				order x2groupi gj vname v*
				sort x2groupi gj
				gen cj=_n
				global cjMax = _N
				drop gj vname
				quietly reshape long v, i(x2groupi cj) j(k)
				gen ri = (x2groupi-1)*$k1+k
				keep ri cj v
				order ri cj v
				sort ri cj v
				by ri: gen t1=_n==1
				quietly expand t1 * ${cjMax}
				sort ri cj v
				by ri cj: gen byte t2= t1==1 & _n>1
				sort ri t2
				by ri t2: replace cj=_n if t2==1 
				replace v = 0 if t2==1
				sort ri cj t2
				by ri cj: drop if t2==1 & _N>1
				keep ri cj v
				quietly reshape wide v, i(ri) j(cj)
				drop ri
				compress
				save "${jtmp}/tMatrixBigMat.dta", replace
				//
				use "${jtmp}/tMatrixBigMat.dta", clear
				mkmat *, matrix(BigMat)
				matrix colnames BigMat = $groupnamesAeq
				matrix rownames BigMat = $groupnamesBeq
				capture erase "${jtmp}/tMatrixBigMat.dta"
				} // Bigmat end
			disp "(3) BigMat produced"
			k_datetime
			/*--------------------------------------------------------*/
			//// Matrix: MV MVlong
			$quietlyFlag {
				// create tvarlistB.dta
				{
					use "${jtmp}/tvarlistA.dta", clear
					quietly keep if inlist(vtype,2,3) 
					sort vi
					gen ri=_n
					save "${jtmp}/tvarlistAa.dta", replace
					//
					global matrixList BetaU
					$quietlyFlag SaveMatrixList matrixList "tMatrix"
					//
					use "${jtmp}/tMatrixBetaU.dta", clear
					xpose, clear
					gen ri= _n
					quietly merge 1:1 ri using "${jtmp}/tvarlistAa.dta", keep(3) nogen
					rename v1 b
					gen j=_n
					keep j vi b x2groupi x2groupnm gj vname
					order j vi x2groupi x2groupnm gj vname b
					save "${jtmp}/tvarlistB.dta", replace		
					//
					capture erase "${jtmp}/tvarlistAa.dta"
					capture erase "${jtmp}/tMatrixBetaU.dta"
					}
				// create MVlong
				{
					matrix MVlong = BetaU
					matrix colnames MVlong = __UN__:			
					}
				// create MV
				{
					capture matrix drop MV
					//
					use "${jtmp}/DecomFigureA.dta", clear
					foreach groupi of numlist 1/$groupN {
						{
							preserve
							use "${jtmp}/tvarlistB.dta", clear
							keep b gj x2groupi vname
							quietly keep if x2groupi==`groupi'
							mkmat b, matrix(b)
							matrix b = b'
							local tv = "groupnamesA`groupi'"
							matrix colnames b = $`tv'
							restore
							} // end preserve/restore pair
						matrix score double Hhat`groupi' = b
						matrix vecaccum hhataccum = Hhat`groupi' $x1all
						matrix coef = hhataccum * X1X1Inv
						matrix score double fitted = coef
						gen double resid`groupi' = Hhat`groupi'-fitted
						matrix MV = nullmat(MV), coef
						matrix drop coef hhataccum
						drop fitted
						} // end foreach groupi
					matrix colnames MV = $groupnamesBeq				
					}
				// update on MVlong
				matrix MVlong = MVlong, MV
				global MVlongNames : colfullnames MVlong
				} // end Matrix: MV MVlong
			// 
			unab tvarlist: resid*
			global residString = "`tvarlist'"
			//
			$quietlyFlag save "${jtmp}/DecomFigureA.dta", replace
			disp "(4) MVlong produced"
			k_datetime	
			/*--------------------------------------------------------*/
			} // end block 6.1
	/*--------------------------------------------------------*/	
		// block 6.2: XXbig, XXBigInv
		{
			matrix TT = J($k1+$k2,$groupN * $k1, 0)
			matrix DD = I($groupN) # X1X1 
			matrix XXBig = (XX, TT) \ (TT', DD)
			matrix drop TT DD
			if 1 {
				matrix TT = J($k1+$k2,$groupN * $k1, 0)
				matrix DD = I($groupN) # X1X1Inv 
				matrix XXBigInv = (XXInv, TT) \ (TT', DD)
				matrix drop TT DD		
				}
			else {
				matrix XXBigInv = syminv(XXBig)
				}
			matrix colnames XXBigInv = $MVlongNames
			matrix rownames XXBigInv = $MVlongNames
			disp "(5) XXBigInv produced"
			k_datetime		
			} // block 6.2 end
	/*--------------------------------------------------------*/	
		// block 6.3: CBigAux, AuxRegV
		{
			use "${jtmp}/DecomFigureA.dta", clear
			if "$clusterV"=="" {
				matrix accum residaccum = ehat $residString, deviations nocons
				matrix residaccum = residaccum / (_N -$k1 - $k2)
				matrix SigmaVBeta = residaccum[2...,1]
				matrix EE = SigmaVBeta # ( (I($k1), Gamma ) * XXInv )
				matrix FF = residaccum[2...,2...] # X1X1Inv
				matrix CBigAux = ( residaccum[1,1]*XXInv, EE' ) \ ( EE, FF )
				matrix colnames CBigAux = $MVlongNames
				matrix rownames CBigAux = $MVlongNames	
				matrix drop EE FF
				global matrixList = "$matrixList residaccum SigmaVBeta"
				}
			else {
				matrix CBigAux = XXBigInv
				matrix colnames CBigAux = $MVlongNames
				matrix rownames CBigAux = $MVlongNames	
				$quietlyFlag _robust ehat $residString, v(CBigAux) cluster($clusterV)
				}
			matrix AuxRegV = CBigAux[1+$k1+$k2..., 1+$k1+$k2...]
			disp "(6) VarCov produced"
			k_datetime
			}
	/*--------------------------------------------------------*/
		// block 6.4: __TC
		{
			use "${jtmp}/DecomFigureA.dta", clear
			matrix score double _tc = BetaUX2
			if "clusterV"~="" {
				quietly reg _tc $x1all, cluster($clusterV)
				}
			else {
				quietly reg _tc $x1all 
				}
			matrix __TCb = e(b)
			matrix coleq __TCb = __TC
			matrix _tcregpart = e(V)
			matrix __TCcov = _tcregpart +Gamma*VU2*Gamma'
			matrix coleq __TCcov = __TC
			matrix roweq __TCcov = __TC	
			}
		}
	/*--------------------------------------------------------*/
	// section 7.9: final organization
	{
		$quietlyFlag {
			use "${jtmp}/tvarlistA.dta", clear
			keep vi x2groupi vtype 
			quietly keep if vtype==3
			keep vi x2groupi 
			quietly replace vi = _n
			sort x2groupi vi
			quietly by x2groupi: keep if inlist(_n,1,_N)
			by x2groupi: gen vj = _n
			quietly reshape wide vi, i(x2groupi) j(vj)
			capture confirm new variable vi2, exact
			if !_rc{
				gen vi2=.
				}
			quietly replace vi2 = vi1 if vi2==.
			//
			matrix VU2Part = J($groupN*$k1,$groupN *$k1,0)
			foreach tgi of numlist 1/$groupN {
				gen byte t1 = 1- (x2groupi==`tgi')
				sort t1
				local k2s= vi1[1]
				local k2e= vi2[1]
				matrix TT = BigMat[1+(`tgi'-1)*$k1..`tgi'*$k1,`k2s'..`k2e']
				matrix SS = VU2[`k2s'..`k2e',`k2s'..`k2e']
				matrix RR = TT * SS * TT'
				matrix VU2Part[1+(`tgi'-1)*$k1,1+(`tgi'-1)*$k1] = RR
				drop t1
				}
			}
		disp "(7) VU2Part produced"
		k_datetime
		/*--------------------------------------------------------*/
		$quietlyFlag {
			clear
			quietly svmat Gamma, names(col)
			quietly xpose, varname clear
			unab tvarlist: v*
			disp "*******: `tvarlist'"
			rename _varname vname
			quietly merge 1:1 vname using "${jtmp}/tvarlistA.dta", keep(3) nogen keepusing(gj x2groupi x2groupnm)
			sort x2groupi gj
			foreach tgi of numlist 1/$groupN {
				quietly mkmat `tvarlist' if x2groupi ==`tgi', matrix(tt) rownames(vname) roweq(x2groupnm)
				matrix Gamma`tgi' = (tt)'
				}
			//
			{
				use "${jtmp}/tvarlistA.dta", clear
				quietly keep if vtype==3
				sort x2groupi gj
				quietly by x2groupi: keep if inlist(_n,1,_N)
				by x2groupi: gen vj = _n
				gen SG = "__UN__:"+vname
				gen SH = x2groupnm+":"
				keep x2groupi SG SH vj
				quietly reshape wide SG SH, i(x2groupi) j(vj)
				capture confirm new variable SG2, exact
				if !_rc{
					gen SG2=""
					gen SH2=""
					}
				quietly replace SG2 = SG1 if SG2==""
				quietly replace SH2 = SH1 if SH2==""
				foreach tgi of numlist 1/$groupN {
					global GroupGS`tgi' = SG1[`tgi']
					global GroupGE`tgi' = SG2[`tgi']
					global GroupH`tgi' = SH1[`tgi']
					}
				}
			//
			capture matrix drop MM
			foreach tgi of numlist 1/$groupN {
				capture matrix drop Mi
				foreach thj of numlist 1/$groupN {
					local gis = "GroupGS`tgi'"
					local gie = "GroupGE`tgi'"
					local hjs = "GroupH`thj'"
					local gamma = "Gamma`tgi'"
					matrix AA = `gamma' * CBigAux["$`gis'".."$`gie'", "$`hjs'"]
					matrix Mi = nullmat(Mi),AA
					}
				matrix MM = nullmat(MM) \ Mi
				}
			matrix Covterms = MM
			}
		disp "(8) Covterms produced"
		//
		matrix FullCov = AuxRegV + VU2Part + Covterms + Covterms'
		matrix colnames FullCov = $groupnamesBeq
		matrix rownames FullCov = $groupnamesBeq
		k_datetime
		}
	/*--------------------------------------------------------*/
	// section 7.10: report
	{
		$quietlyFlag {
		matrix TT = FullCov \ MV
		//
		use "${jtmp}/DecomFigureA.dta", clear
		global nobs = _N
		//
		clear
		quietly svmat FullCov, names(v)
		gen ri = _n
		quietly reshape long v, i(ri) j(j)
		gen tri = mod(ri-1,$k1)+1
		gen tj = mod(j-1,$k1)+1
		sort tri tj ri j
		quietly drop if tri == $k1 | tj==$k1
		sort ri j
		quietly by ri: replace j=_n
		sort j ri
		quietly by j: replace ri=_n
		keep ri j v
		quietly reshape wide v, i(ri) j(j)
		drop ri
		quietly mkmat v*, matrix(GCov)
		//matrix list GCov
		local tgn: word count $x2deltaVarlist
		if _N > `tgn' {
			matrix colnames GCov = $x2deltaVarlist Residual
			matrix rownames GCov = $x2deltaVarlist Residual
			}
		else {
			matrix colnames GCov = $x2deltaVarlist 
			matrix rownames GCov = $x2deltaVarlist 
			}
		//
		clear
		svmat MV, names(v)
		xpose, clear
		gen ri=_n
		gen tri = mod(ri-1,$k1)+1
		quietly drop if tri == $k1
		quietly replace ri=_n
		mkmat v1, matrix(Delta)
		matrix Delta = Delta'
		local tgn: word count $x2deltaVarlist
		if _N > `tgn' {
			matrix colnames Delta = $x2deltaVarlist Residual
			}
		else {
			matrix colnames Delta = $x2deltaVarlist 
			}
		//matrix list Delta	
		matrix R = I($groupN) \ J(1,$groupN, 1)
		matrix DeltaB = Delta*R'
		matrix GCovB = R*GCov*R'
		local tgn: word count $x2deltaVarlist
		if _N > `tgn' {
			matrix colnames DeltaB = $x2deltaVarlist Residual _TC
			}
		else {
			matrix colnames DeltaB = $x2deltaVarlist _TC
			}
		local tgn: word count $x2deltaVarlist
		if _N > `tgn' {
			matrix colnames GCovB = $x2deltaVarlist Residual _TC
			matrix rownames GCovB = $x2deltaVarlist Residual _TC
			}
		else {
			matrix colnames GCovB = $x2deltaVarlist _TC
			matrix rownames GCovB = $x2deltaVarlist _TC
			}		
		//
		matrix DeltaBB = DeltaB
		matrix GCovBB = GCovB
		}
		//
		disp "(9) report...."
		ereturn post DeltaBB GCovBB, obs(${nobs}) depname($depvar) 
		ereturn display
		macro list depvar x1all 
		capture macro list ClusterV
		k_datetime
		}
	/*--------------------------------------------------------*/
	// section 7.11: prepare for coefficient plot
	{
		$quietlyFlag {
			matrix t = DeltaB \ vecdiag(GCovB)
			clear
			quietly svmat t, names(col)
			quietly xpose, varname clear
			rename _varname vname
			order vname v1 v2
			rename v1 b
			rename v2 v
			gen se = sqrt(v)
			gen range = se*invnormal(0.95)
			gen lb = b-range
			gen ub = b+range
			drop if vname=="_TC"
			mkmat b, matrix(bb) rownames(vname)
			mkmat lb ub, matrix(ci) rownames(vname)
			matrix figB = bb'
			matrix figCI = ci'
			matrix drop bb ci t
			}
		}
	/*--------------------------------------------------------*/
	global matrixList DeltaB GCovB figB figCI
	macro drop matrixListSave 
	foreach tv of global matrixList {
		matrix rename `tv' `tv'`tSeq'
		global matrixListSave = "$matrixListSave `tv'`tSeq'"
		}
	quietly SaveMatrixList matrixListSave "tMatrix"
	macro drop matrixListSave matrixList
	//
	use "${jtmp}/tvarlistA.dta", clear
	save "${jtmp}/tvarlistA`tSeq'.dta", replace
	//
	capture erase "${jtmp}/tvarlistA.dta"
	capture erase "${jtmp}/tvarlistB.dta"
	//
	disp "(10) done"
	k_datetime
	end
	}
/*--------------------------------------------------------*/