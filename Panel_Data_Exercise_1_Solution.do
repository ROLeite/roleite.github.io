*==============================================================================
* Panel Data Exercise 1 -- Solution do-file
* Pooled OLS, Random Effects (FGLS) and Fixed Effects
* Data: Cornwell and Rupert Returns to Schooling Data (595 individuals, T = 7)
*==============================================================================
* Notation: the exercise writes eps_it = c_i + u_it.
*   Stata's e(sigma_u) = sd of the individual effect c_i      (exercise: sigma_c)
*   Stata's e(sigma_e) = sd of the idiosyncratic error u_it  (exercise: sigma_u)
*==============================================================================

clear all
set more off
version 15

*------------------------------------------------------------------------------
* 0. Data
*------------------------------------------------------------------------------
* Change the path/file name to the copy downloaded from the course website.
local datafile "cornwell_rupert.csv"

import delimited using "`datafile'", clear case(lower)
* (if you have a .dta file instead: use "cornwell_rupert.dta", clear
*  followed by: rename *, lower)

* Individual and time identifiers (data are sorted by individual, 7 rows each)
capture confirm variable id
if _rc gen id = ceil(_n/7)
capture drop t
bysort id: gen t = _n

xtset id t
xtdescribe

global X  occ smsa ms fem ed exp      // full regressor list
global XV occ smsa ms exp             // time-varying regressors

summarize lwage $X

*------------------------------------------------------------------------------
* PART A -- Pooled OLS
*------------------------------------------------------------------------------
* Q1. POLS with conventional and clustered standard errors
regress lwage $X
estimates store pols
scalar SSR_pols = e(rss)

regress lwage $X, vce(cluster id)
estimates store pols_cl

estimates table pols pols_cl, b(%9.4f) se(%9.4f) stats(N r2)

* Q2. FEM coefficient: exact percentage gap and t tests
quietly regress lwage $X
display "FEM coefficient = " %7.4f _b[fem]
display "Exact % wage gap = " %6.2f 100*(exp(_b[fem]) - 1) "%"
test fem = 0                                  // conventional SE

quietly regress lwage $X, vce(cluster id)
test fem = 0                                  // clustered SE

*------------------------------------------------------------------------------
* PART B -- Should we use GLS?
*------------------------------------------------------------------------------
* Q4a. Correlation of POLS residuals within individual (t and t-1)
quietly regress lwage $X
predict e_pols, residuals
correlate e_pols L.e_pols

* Q4b. Breusch-Pagan LM test, by hand ...
bysort id: egen sum_e = total(e_pols)
gen e2 = e_pols^2
quietly summarize e2
scalar SSE = r(sum)
egen tag = tag(id)
gen sum_e2 = sum_e^2 if tag
quietly summarize sum_e2
scalar NT = _N
scalar T  = 7
scalar LM = NT/(2*(T-1)) * (r(sum)/SSE - 1)^2
display "Breusch-Pagan LM = " %9.2f LM "   p-value = " %6.4f chi2tail(1, LM)
drop sum_e e2 sum_e2

* ... and with Stata's command (after RE)
quietly xtreg lwage $X, re
xttest0

*------------------------------------------------------------------------------
* PART C -- Random effects (FGLS)
*------------------------------------------------------------------------------
* Q5. RE estimates, variance components, rho and theta
xtreg lwage $X, re theta
estimates store re
scalar sig2_c = e(sigma_u)^2        // individual effect (exercise: sigma_c^2)
scalar sig2_u = e(sigma_e)^2        // idiosyncratic     (exercise: sigma_u^2)
scalar theta  = e(theta)
display "sigma_c^2 = " %7.4f sig2_c "   sigma_u^2 = " %7.4f sig2_u
display "rho       = " %7.4f e(rho)  "   theta     = " %7.4f theta
display "theta by formula = " %7.4f 1 - sqrt(sig2_u/(sig2_u + T*sig2_c))

* Q6. RE as OLS on quasi-demeaned data
foreach v of varlist lwage $X {
    bysort id: egen m_`v' = mean(`v')
    gen qd_`v' = `v' - theta*m_`v'
}
gen qd_const = 1 - theta
regress qd_lwage qd_occ qd_smsa qd_ms qd_fem qd_ed qd_exp qd_const, noconstant
* -> coefficients identical to xtreg, re (standard errors differ slightly by dof)

* Q7. Compare FEM and ED: POLS vs RE
estimates table pols re, keep(fem ed) b(%9.4f) se(%9.4f)

*------------------------------------------------------------------------------
* PART D -- Fixed effects
*------------------------------------------------------------------------------
* Q8. Within estimator and LSDV -- identical slopes; FEM and ED are dropped
xtreg lwage $X, fe
estimates store fe

areg lwage $X, absorb(id)            // LSDV (individual dummies absorbed)
estimates store lsdv
scalar SSR_lsdv = e(rss)
* equivalent, explicit dummies (slower):  regress lwage $X i.id

estimates table fe lsdv, keep($XV) b(%9.4f) se(%9.4f)

* Q9. F test of no individual effects, by hand
quietly xtreg lwage $XV, fe
scalar N  = e(N_g)                   // number of individuals (595)
scalar K  = 4                        // time-varying slopes in FE
scalar df2 = NT - N - K              // residual df of LSDV (3,566)

* (a) Restricted model = POLS from Q1 (includes FEM and ED)
*     LSDV nests FEM and ED in the dummies -> N - 3 restrictions
scalar q  = N - 3
scalar F  = ((SSR_pols - SSR_lsdv)/q) / (SSR_lsdv/df2)
display "F = " %8.3f F "   df = (" q ", " df2 ")   p-value = " %6.4f Ftail(q, df2, F)

* (b) Restricted model = POLS with time-varying regressors only
*     -> N - 1 restrictions; this is the test printed by -xtreg, fe-
quietly regress lwage $XV
scalar SSR_tv = e(rss)
scalar q  = N - 1
scalar F  = ((SSR_tv - SSR_lsdv)/q) / (SSR_lsdv/df2)
display "F = " %8.3f F "   df = (" q ", " df2 ")   p-value = " %6.4f Ftail(q, df2, F)
display "5% critical value = " %6.3f invFtail(q, df2, 0.05)

* Q10. FE with year dummies: EXP is dropped (collinear with the time trend)
xtreg lwage $X i.t, fe
* Check: exp - t is constant within individual
gen exp_minus_t = exp - t
bysort id: egen sd_check = sd(exp_minus_t)
summarize sd_check                    // all zeros
drop exp_minus_t sd_check

* Q11. Between estimator and comparison
xtreg lwage $X, be
estimates store be

estimates table pols be re fe, keep($XV fem ed) b(%9.4f) se(%9.4f)

*------------------------------------------------------------------------------
* PART E -- Fixed or random effects?
*------------------------------------------------------------------------------
* Q12. Classic Hausman test
hausman fe re
hausman fe re, sigmamore              // uses a common sigma

* Is V_FE - V_RE positive definite? (compare the 4 time-varying coefficients)
quietly xtreg lwage $XV, fe
matrix bF = e(b)
matrix VF = e(V)
quietly xtreg lwage $XV fem ed, re    // same order: occ smsa ms exp first
matrix bR = e(b)
matrix VR = e(V)
mata:
    bF = st_matrix("bF")[1, 1..4]
    bR = st_matrix("bR")[1, 1..4]
    D  = st_matrix("VF")[1..4, 1..4] - st_matrix("VR")[1..4, 1..4]
    printf("Eigenvalues of V_FE - V_RE:\n")
    symeigenvalues(D)
    d  = (bF - bR)'
    H  = d' * invsym(D) * d
    printf("Hausman statistic (by hand) = %9.2f\n", H)
end
* -> negative eigenvalues: the classic statistic is unreliable; use Q13

* Q13. Mundlak / correlated random effects (robust Hausman)
* (individual means m_occ m_smsa m_ms m_exp were created in Q6)
xtreg lwage $X m_occ m_smsa m_ms m_exp, re vce(cluster id)
estimates store cre

* (a) Coefficients on time-varying regressors equal FE
estimates table fe cre, keep($XV) b(%9.4f)

* (b) Robust Hausman: are the individual means jointly zero? chi2(4)
*     Rejecting means c_i is correlated with the regressors -> RE inconsistent
estimates restore cre
test m_occ m_smsa m_ms m_exp

* (c) Preferred model: CRE/FE. CRE still estimates the FEM gap
display "CRE FEM coefficient = " %7.4f _b[fem]
display "CRE FEM % gap       = " %6.2f 100*(exp(_b[fem]) - 1) "%"

*------------------------------------------------------------------------------
* PART F -- Different wage equations for men and women
*------------------------------------------------------------------------------
* Q14. Separate RE regressions
preserve
    keep if tag
    tab fem                           // number of individuals by sex
restore
tab fem                               // number of observations by sex

xtreg lwage occ smsa ms ed exp if fem == 0, re
estimates store re_men
xtreg lwage occ smsa ms ed exp if fem == 1, re
estimates store re_women
estimates table re_men re_women, b(%9.4f) se(%9.4f)

* Q15a. Pooled RE with full interactions, robust Wald (Chow-type) test
foreach v in occ smsa ms ed exp {
    gen f_`v' = fem*`v'
}
xtreg lwage occ smsa ms ed exp fem f_occ f_smsa f_ms f_ed f_exp, ///
    re vce(cluster id)
test fem f_occ f_smsa f_ms f_ed f_exp

* Q15b. FE with interactions of FEM and time-varying regressors
*       (FEM and FEM*ED are time invariant -> not identified)
xtreg lwage occ smsa ms exp f_occ f_smsa f_ms f_exp, fe vce(cluster id)
test f_occ f_smsa f_ms f_exp

*==============================================================================
* End of file
*==============================================================================
