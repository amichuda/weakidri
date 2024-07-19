*! v1.1.1 - Small change to how results are opened after running inference procedure
*! v1.1.0 - Added ability to include endogenous covariates
*! v1.0.0 - Aleksandr Michuda, am2497@cornell.edu

program define weakidri, rclass
    version 16.0
    syntax varname [if], hhid(varname) Hybrid(varname) MIn(real) MAx(real) INCrement(real) [STORE_results(name) progress(integer 1) path(string) TEST_type(string) BASEgroup(integer 2) type_one(real 0.05) ADD_post(string) controls(varlist fv) ENDOgenous(varname) time(varname) endo_factor]
    /*
        Program that carries out a brute force search for phi for an outcome `varname', with household ID `hhid', and
        hybrid variable `hybrid', and optionally, controls `controls'.

        The program creates a grid from `min' to `max' with an increment of `increment'. Optionally you can set the
        Type-I error with `type_one'.

        This program saves p-values to a postfile and you can add another variable or local, with `add_post'.

        There are three types of tests that are implemented (set with `test_type'):
        1. "joint": This tests all restrictions as a joint test. This is the default setting.
        2. "separate": This tests all restrictions separately and p-values are saved separately. It this is chosen for T=2,
            it switches to "joint" as those are equivalent in the T=2 case.
        3. "base" : This test separate restrictions against a base group `base_group', set by the number of the trajectory.
            (run tab trajectory, to see which number to use)

        Store estimates by adding a name to `store_results'

        if `progress' is chosen, will omit output from the test command and show dots for each test that is run.

        Optionally save to `path'

        Args:
            varname: The dependent variable to be used in the regression
            hhid(varname): The household id. Used for clustering SE for auxiliary regression
            hybrid(varname): the binary variable of adoption
            min(real): The left endpoint of the grid
            max(real): The right endpoint of the grid
            increment(real): How finely to break up the grid
            store_results(name, optional): use store the auxiliary regression estimates with `name' for later use
            progress(integer, optional): Whether to show the output of each test that is run or show dots on the screen. Default to true.
            path(string, optional): where to save the results of the simulation
            test_type(string, optional): The type of test to run, described above.
            basegroup(integer, optional): the basegroup to use with test_type "base"
            type_one(real, optional): The Type-I error of the test. Defaults to 0.05
            add_post(string, optional): Optionally save another parameter to the simulation dataset. Only for use with `path'
            controls(varlist, optional): Optionally include controls into the regression
            endogenous(varname, optional): Optionally include an endogenous covariate to co-determine theta along with adoption.
            endo_factor (optionally on, optional): Optionally specify that the endogenous variable is a factor variable so that the regression will use `i` syntax. 
                Otherwise, a continuous variable is assumed.

        Returns:
            r(min_phi`i'`j'): The minimum phi that failed to reject the test; the left side of the confidence interval for `phi', for switching trajectories `i' and `j'
            r(max_phi`i'`j'): The maximum phi that failed to reject the test; the right side of the confidence interval for `phi', for switching trajectories `i' and `j'
            r(min_phi_joint): The minimum phi that failed to reject the test of the joint test of all switching trajectories
            r(max_phi_joint): The maximum phi that failed to reject the test of the joint test of all switching trajectories
            r(table): The table of regression results


    */

    * Get post file for saving p-values
	tempname memhold

    if "`path'" != "" {
        local results_data "`path'"
    }
    else {
        tempfile results_phi
        local results_data "`results_phi'"
    }

    if "`test_type'" == "base" {
        local memhold_vars `"phi p_val`basegroup'3 p_val`basegroup'4 p_val`basegroup'5 p_val`basegroup'6 p_val`basegroup'7" p_val_joint'
        if "`endogenous'" != "" {
            local memhold_vars `"`memhold_vars' p_val_endo_1`basegroup'3 p_val_endo_1`basegroup'4 p_val_endo_1`basegroup'5 p_val_endo_1`basegroup'6 p_val_endo_1`basegroup'7"'
            local memhold_vars `"`memhold_vars' p_val_endo_2`basegroup'3 p_val_endo_2`basegroup'4 p_val_endo_2`basegroup'5 p_val_endo_2`basegroup'6 p_val_endo_2`basegroup'7"'

        }
    }
    else if "`test_type'" == "joint" {
        local memhold_vars `"phi p_val_joint"'
    }
    else {
        local memhold_vars `"phi p_val32 p_val43 p_val54 p_val65 p_val76"'
        if "`endogenous'" != "" {
            local memhold_vars `"`memhold_vars' p_val_endo_1_32 p_val_endo_1_43 p_val_endo_1_54 p_val_endo_1_65 p_val_endo_1_76"'
            local memhold_vars `"`memhold_vars' p_val_endo_2_32 p_val_endo_2_43 p_val_endo_2_54 p_val_endo_2_65 p_val_endo_2_76"'
        }
    }

    postfile `memhold' `memhold_vars' str20(test_type add) using "`results_data'", replace


    local          never 1
    tab             trajectory

    local         always `r(r)'
    local          lastswitcher = `always'-1

    numlist         "2(1)`lastswitcher'"
    local          switchers `r(numlist)'

    numlist         "0(1)`lastswitcher'"
    local          noalways `r(numlist)'

    * expand factor variables if needed
    local fvops = "`s(fvops)'" == "true" | _caller() >= 11 & "`controls'" != ""
    if `fvops' {
        fvexpand `controls'
        local controls =  r(varlist)
    }
    else {
        local controls
    }

    if "`endo_factor'" == "endo_factor" {
        local ef  i
    }
    else {
        local ef  c
    }

    * If endo is included, widen variable repeated by household
    if "`endogenous'" != "" {
        if "`time'" == "" {
            di as err "if endogenous is specified, a time variable must also be specified."
            exit 102
        }
        preserve
            keep `hhid' `time' `endogenous'
            reshape wide `endogenous', i(`hhid') j(`time')
            * save new variable names
            local endo_wide_vars: char _dta[ReS_Xij_wide1]
            tempfile endo_
            save `endo_'
        restore

        merge m:1 `hhid' using `endo_'
        drop _merge
    }

    if "`endogenous'" == "" {
         reg    `varlist' i(`noalways').trajectory                ///
                i(`switchers' `always').trajectory#1.`hybrid' 	         ///
                `controls' ///
                `if', vce(cluster `hhid') nocons

    }

    else {
        * create interactions with new endogenous variable
        local endo_interactions 

        foreach endo_var of local endo_wide_vars {
            local endo_interactions `endo_interactions' i(`noalways').trajectory#`ef'.`endo_var' i(`switchers'  `always').trajectory#1.`hybrid'#`ef'.`endo_var'
        }       

        reg    `varlist' i(`noalways').trajectory i(`switchers'  `always').trajectory#1.`hybrid' ///
        `controls' `endogenous' `endo_interactions' ///
        `if', vce(cluster `hhid') nocons
    }

    if "`store_results'" != "" {
        estimates store `store_results'
    }


    qui ereturn list
    return add

	/*
	We test for a plethora of phi values
	*/
    if `progress' == 1 {
        _dots 0, title("Running Weak-ID Robust Inference Test for hhid: `hhid', controls: `controls' and endogenous variables: `endogenous', with min: `min', max: `max' and increment: `increment'")
    }


	forvalues phi_potential = `=`min''(`=`increment'')`=`max'' {

        if `progress' == 1 {
            local i = `i' + 1
            _dots `i' 0
            local qui quietly
        }

        qui tab trajectory
        local test_string 
        if r(r) == 4 & "`endogenous'" == "" {
            local test_type = "joint"
            * If T=2, with no endogenous covariates, joint and separate are the same
            * No longer the case, but keeping here for posterity for
        }
        if "`test_type'" == "joint" {
            forvalues i = `=`r(r)'-1'(-1)3 {
                local test_string ((`i'.trajectory#1.`hybrid'- `=`i'-1'.trajectory#1.`hybrid') = `phi_potential'*(`i'.trajectory-`=`i' -1'.trajectory)) `test_string'
            
                if "`endogenous'" != "" {
                    * add endogenous covariates to test_string
                    foreach endo_var of local endo_wide_vars {
                        local test_string ((`i'.trajectory#1.`hybrid'#`ef'.`endo_var' - `=`i'-1'.trajectory#1.`hybrid'#`ef'.`endo_var') = `phi_potential'*(`i'.trajectory#`ef'.`endo_var' -`=`i' -1'.trajectory#`ef'.`endo_var')) `test_string'
                    }
                }
            }
            `qui' test `test_string'
            qui return list

            loc p_v = r(p)
            post `memhold' (`phi_potential') (`p_v') ("`test_type'") ("`add_post'") 

        }
        else if "`test_type'" == "separate" {

            forvalues i = `=`r(r)'-1'(-1)3 {
                `qui' test ((`i'.trajectory#1.`hybrid'-`=`i'-1'.trajectory#1.`hybrid') = `phi_potential'*(`i'.trajectory-`=`i' -1'.trajectory))
                scalar p_v_`i' = r(p)
            }

            local memhold_save `"(p_v_3) (p_v_4) (p_v_5) (p_v_6) (p_v_7)"'

            if "`endogenous'" != "" {
                * add endogenous covariates to test_string
                foreach endo_var of local endo_wide_vars {
                    forvalues i = `=`r(r)'-1'(-1)3 {
                        `qui' test ((`i'.trajectory#1.`hybrid'#`ef'.`endo_var' - `=`i'-1'.trajectory#1.`hybrid'#`ef'.`endo_var') = `phi_potential'*(`i'.trajectory#`ef'.`endo_var' -`=`i' -1'.trajectory#`ef'.`endo_var')) `test_string'
                        scalar p_v_`i'_`endo_var' = r(p)
                    }
                    local memhold_save `"`memhold_save' (p_v_3_`endo_var') (p_v_4_`endo_var') (p_v_5_`endo_var') (p_v_6_`endo_var') (p_v_7_`endo_var')"' 
                }

            }
            post `memhold' (`phi_potential') `memhold_save' ("`test_type'") ("`add_post'")
        }

        else if "`test_type'" == "base" {
            if "basegroup" == "" {
                di "Base group not chosen. Choosing trajectory 1"
                local basegroup 1
            }

            if "`endogenous'" != "" {
                di in error "Using a trajectory as base with endogenous variables not implemented."
                exit 187
            }

            local test_string 

            forvalues i = 2(1)`=`r(r)'-1' {
                if `i' == `basegroup' {
                    continue
                }
                local test_string ((`i'.trajectory#1.`hybrid'-`basegroup'.trajectory#1.`hybrid') = `phi_potential'*(`i'.trajectory-`basegroup'.trajectory)) `test_string'
                `qui' test ((`i'.trajectory#1.`hybrid'-`basegroup'.trajectory#1.`hybrid') = `phi_potential'*(`i'.trajectory-`basegroup'.trajectory))
                scalar p_v_`i' = r(p)
            }

            `qui' test `test_string'
            qui return list
            scalar p_val_joint = r(p)
            post `memhold' (`phi_potential') (p_v_3) (p_v_4) (p_v_5) (p_v_6) (p_v_7) (p_val_joint) ("`test_type'") ("`add_post'")

        }

	}
    tab trajectory

    postclose `memhold'

    preserve
        use "`results_data'", clear
        destring *, replace

        if test_type == "joint" {
            gen not_sig = (p_val_joint>=`type_one')
            su phi if not_sig == 1

            gen min_phi = round(r(min),0.001)
            gen max_phi = round(r(max),0.001)

            loc min_phi = round(r(min),0.001)
            loc max_phi = round(r(max),0.001)
            loc mean_phi = r(mean)

            return local min_phi = `min_phi'
            return local max_phi = `max_phi'
        }

        if test_type == "separate" {
            gen not_sig32 = (p_val32>=`type_one')
            su phi if not_sig32 == 1

            gen min_phi32 = round(r(min),0.001)
            gen max_phi32 = round(r(max),0.001)

            loc min_phi32 = round(r(min),0.001)
            loc max_phi32 = round(r(max),0.001)
            loc mean_phi32 = r(mean)

            return local min_phi32 = `min_phi32'
            return local max_phi32 = `max_phi32'

            forval i=4/7 {
                gen not_sig`i'`=`i'-1' = (p_val`i'`=`i'-1' >= `type_one')
                su phi if not_sig`i'`=`i'-1' ==1
                loc min_phi`i'`=`i'-1' = round(r(min), 0.001)
                loc max_phi`i'`=`i'-1' = round(r(max), 0.001)

                * Also create variables for those terms
                gen min_phi`i'`=`i'-1' = round(r(min), 0.001)
                gen max_phi`i'`=`i'-1' = round(r(max), 0.001)

                return local min_phi`i'`=`i'-1' = `min_phi`i'`=`i'-1''
                return local max_phi`i'`=`i'-1' = `max_phi`i'`=`i'-1''
            }
        }

        else if test_type == "base" {
            forval i=3/7 {
                gen not_sig2`i' = (p_val2`i' >= `type_one')
                su phi if not_sig2`i' ==1
                loc min_phi2`i' = round(r(min), 0.001)`'
                loc max_phi2`i' = round(r(max), 0.001)

                * Also create variables for those terms
                gen min_phi2`i' = round(r(min), 0.001)
                gen max_phi2`i' = round(r(max), 0.001)


                return local min_phi2`i' = `min_phi2`i''
                return local max_phi2`i' = `max_phi2`i''
            }
        }

        if "`path'" != "" {
            save "`path'", replace
        }

    restore


end
