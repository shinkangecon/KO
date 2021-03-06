%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'kang_linear';
M_.dynare_version = '5.0';
oo_.dynare_version = '5.0';
options_.dynare_version = '5.0';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(10,1);
M_.exo_names_tex = cell(10,1);
M_.exo_names_long = cell(10,1);
M_.exo_names(1) = {'e_v'};
M_.exo_names_tex(1) = {'e\_v'};
M_.exo_names_long(1) = {'e_v'};
M_.exo_names(2) = {'e_xi'};
M_.exo_names_tex(2) = {'e\_xi'};
M_.exo_names_long(2) = {'e_xi'};
M_.exo_names(3) = {'e_phi'};
M_.exo_names_tex(3) = {'e\_phi'};
M_.exo_names_long(3) = {'e_phi'};
M_.exo_names(4) = {'e_psi'};
M_.exo_names_tex(4) = {'e\_psi'};
M_.exo_names_long(4) = {'e_psi'};
M_.exo_names(5) = {'e_z'};
M_.exo_names_tex(5) = {'e\_z'};
M_.exo_names_long(5) = {'e_z'};
M_.exo_names(6) = {'e_R'};
M_.exo_names_tex(6) = {'e\_R'};
M_.exo_names_long(6) = {'e_R'};
M_.exo_names(7) = {'e_g'};
M_.exo_names_tex(7) = {'e\_g'};
M_.exo_names_long(7) = {'e_g'};
M_.exo_names(8) = {'e_p'};
M_.exo_names_tex(8) = {'e\_p'};
M_.exo_names_long(8) = {'e_p'};
M_.exo_names(9) = {'e_w'};
M_.exo_names_tex(9) = {'e\_w'};
M_.exo_names_long(9) = {'e_w'};
M_.exo_names(10) = {'e_n'};
M_.exo_names_tex(10) = {'e\_n'};
M_.exo_names_long(10) = {'e_n'};
M_.endo_names = cell(26,1);
M_.endo_names_tex = cell(26,1);
M_.endo_names_long = cell(26,1);
M_.endo_names(1) = {'c'};
M_.endo_names_tex(1) = {'c'};
M_.endo_names_long(1) = {'c'};
M_.endo_names(2) = {'h'};
M_.endo_names_tex(2) = {'h'};
M_.endo_names_long(2) = {'h'};
M_.endo_names(3) = {'l'};
M_.endo_names_tex(3) = {'l'};
M_.endo_names_long(3) = {'l'};
M_.endo_names(4) = {'w'};
M_.endo_names_tex(4) = {'w'};
M_.endo_names_long(4) = {'w'};
M_.endo_names(5) = {'R'};
M_.endo_names_tex(5) = {'R'};
M_.endo_names_long(5) = {'R'};
M_.endo_names(6) = {'pi'};
M_.endo_names_tex(6) = {'pi'};
M_.endo_names_long(6) = {'pi'};
M_.endo_names(7) = {'q'};
M_.endo_names_tex(7) = {'q'};
M_.endo_names_long(7) = {'q'};
M_.endo_names(8) = {'r_k'};
M_.endo_names_tex(8) = {'r\_k'};
M_.endo_names_long(8) = {'r_k'};
M_.endo_names(9) = {'k'};
M_.endo_names_tex(9) = {'k'};
M_.endo_names_long(9) = {'k'};
M_.endo_names(10) = {'x'};
M_.endo_names_tex(10) = {'x'};
M_.endo_names_long(10) = {'x'};
M_.endo_names(11) = {'y'};
M_.endo_names_tex(11) = {'y'};
M_.endo_names_long(11) = {'y'};
M_.endo_names(12) = {'N'};
M_.endo_names_tex(12) = {'N'};
M_.endo_names_long(12) = {'N'};
M_.endo_names(13) = {'u'};
M_.endo_names_tex(13) = {'u'};
M_.endo_names_long(13) = {'u'};
M_.endo_names(14) = {'Omega'};
M_.endo_names_tex(14) = {'Omega'};
M_.endo_names_long(14) = {'Omega'};
M_.endo_names(15) = {'theta_p'};
M_.endo_names_tex(15) = {'theta\_p'};
M_.endo_names_long(15) = {'theta_p'};
M_.endo_names(16) = {'pi_w'};
M_.endo_names_tex(16) = {'pi\_w'};
M_.endo_names_long(16) = {'pi_w'};
M_.endo_names(17) = {'theta_w'};
M_.endo_names_tex(17) = {'theta\_w'};
M_.endo_names_long(17) = {'theta_w'};
M_.endo_names(18) = {'v'};
M_.endo_names_tex(18) = {'v'};
M_.endo_names_long(18) = {'v'};
M_.endo_names(19) = {'xi'};
M_.endo_names_tex(19) = {'xi'};
M_.endo_names_long(19) = {'xi'};
M_.endo_names(20) = {'phi'};
M_.endo_names_tex(20) = {'phi'};
M_.endo_names_long(20) = {'phi'};
M_.endo_names(21) = {'psi'};
M_.endo_names_tex(21) = {'psi'};
M_.endo_names_long(21) = {'psi'};
M_.endo_names(22) = {'z'};
M_.endo_names_tex(22) = {'z'};
M_.endo_names_long(22) = {'z'};
M_.endo_names(23) = {'g'};
M_.endo_names_tex(23) = {'g'};
M_.endo_names_long(23) = {'g'};
M_.endo_names(24) = {'eps_p'};
M_.endo_names_tex(24) = {'eps\_p'};
M_.endo_names_long(24) = {'eps_p'};
M_.endo_names(25) = {'eps_w'};
M_.endo_names_tex(25) = {'eps\_w'};
M_.endo_names_long(25) = {'eps_w'};
M_.endo_names(26) = {'eps_n'};
M_.endo_names_tex(26) = {'eps\_n'};
M_.endo_names_long(26) = {'eps_n'};
M_.endo_partitions = struct();
M_.param_names = cell(71,1);
M_.param_names_tex = cell(71,1);
M_.param_names_long = cell(71,1);
M_.param_names(1) = {'gamma'};
M_.param_names_tex(1) = {'gamma'};
M_.param_names_long(1) = {'gamma'};
M_.param_names(2) = {'zeta'};
M_.param_names_tex(2) = {'zeta'};
M_.param_names_long(2) = {'zeta'};
M_.param_names(3) = {'vartheta'};
M_.param_names_tex(3) = {'vartheta'};
M_.param_names_long(3) = {'vartheta'};
M_.param_names(4) = {'beta'};
M_.param_names_tex(4) = {'beta'};
M_.param_names_long(4) = {'beta'};
M_.param_names(5) = {'delta'};
M_.param_names_tex(5) = {'delta'};
M_.param_names_long(5) = {'delta'};
M_.param_names(6) = {'kappa_x'};
M_.param_names_tex(6) = {'kappa\_x'};
M_.param_names_long(6) = {'kappa_x'};
M_.param_names(7) = {'rho'};
M_.param_names_tex(7) = {'rho'};
M_.param_names_long(7) = {'rho'};
M_.param_names(8) = {'alpha'};
M_.param_names_tex(8) = {'alpha'};
M_.param_names_long(8) = {'alpha'};
M_.param_names(9) = {'f'};
M_.param_names_tex(9) = {'f'};
M_.param_names_long(9) = {'f'};
M_.param_names(10) = {'rho_N'};
M_.param_names_tex(10) = {'rho\_N'};
M_.param_names_long(10) = {'rho_N'};
M_.param_names(11) = {'kappa_u'};
M_.param_names_tex(11) = {'kappa\_u'};
M_.param_names_long(11) = {'kappa_u'};
M_.param_names(12) = {'omega'};
M_.param_names_tex(12) = {'omega'};
M_.param_names_long(12) = {'omega'};
M_.param_names(13) = {'s_m'};
M_.param_names_tex(13) = {'s\_m'};
M_.param_names_long(13) = {'s_m'};
M_.param_names(14) = {'kappa_p'};
M_.param_names_tex(14) = {'kappa\_p'};
M_.param_names_long(14) = {'kappa_p'};
M_.param_names(15) = {'chi_y'};
M_.param_names_tex(15) = {'chi\_y'};
M_.param_names_long(15) = {'chi_y'};
M_.param_names(16) = {'eta_y'};
M_.param_names_tex(16) = {'eta\_y'};
M_.param_names_long(16) = {'eta_y'};
M_.param_names(17) = {'kappa_w'};
M_.param_names_tex(17) = {'kappa\_w'};
M_.param_names_long(17) = {'kappa_w'};
M_.param_names(18) = {'chi_l'};
M_.param_names_tex(18) = {'chi\_l'};
M_.param_names_long(18) = {'chi_l'};
M_.param_names(19) = {'eta_l'};
M_.param_names_tex(19) = {'eta\_l'};
M_.param_names_long(19) = {'eta_l'};
M_.param_names(20) = {'gamma_l'};
M_.param_names_tex(20) = {'gamma\_l'};
M_.param_names_long(20) = {'gamma_l'};
M_.param_names(21) = {'gamma_p'};
M_.param_names_tex(21) = {'gamma\_p'};
M_.param_names_long(21) = {'gamma_p'};
M_.param_names(22) = {'gamma_w'};
M_.param_names_tex(22) = {'gamma\_w'};
M_.param_names_long(22) = {'gamma_w'};
M_.param_names(23) = {'rho_R'};
M_.param_names_tex(23) = {'rho\_R'};
M_.param_names_long(23) = {'rho_R'};
M_.param_names(24) = {'alpha_pi'};
M_.param_names_tex(24) = {'alpha\_pi'};
M_.param_names_long(24) = {'alpha_pi'};
M_.param_names(25) = {'alpha_y'};
M_.param_names_tex(25) = {'alpha\_y'};
M_.param_names_long(25) = {'alpha_y'};
M_.param_names(26) = {'sigma_R'};
M_.param_names_tex(26) = {'sigma\_R'};
M_.param_names_long(26) = {'sigma_R'};
M_.param_names(27) = {'rho_v'};
M_.param_names_tex(27) = {'rho\_v'};
M_.param_names_long(27) = {'rho_v'};
M_.param_names(28) = {'sigma_v'};
M_.param_names_tex(28) = {'sigma\_v'};
M_.param_names_long(28) = {'sigma_v'};
M_.param_names(29) = {'rho_xi'};
M_.param_names_tex(29) = {'rho\_xi'};
M_.param_names_long(29) = {'rho_xi'};
M_.param_names(30) = {'sigma_xi'};
M_.param_names_tex(30) = {'sigma\_xi'};
M_.param_names_long(30) = {'sigma_xi'};
M_.param_names(31) = {'rho_phi'};
M_.param_names_tex(31) = {'rho\_phi'};
M_.param_names_long(31) = {'rho_phi'};
M_.param_names(32) = {'sigma_phi'};
M_.param_names_tex(32) = {'sigma\_phi'};
M_.param_names_long(32) = {'sigma_phi'};
M_.param_names(33) = {'rho_psi'};
M_.param_names_tex(33) = {'rho\_psi'};
M_.param_names_long(33) = {'rho_psi'};
M_.param_names(34) = {'sigma_psi'};
M_.param_names_tex(34) = {'sigma\_psi'};
M_.param_names_long(34) = {'sigma_psi'};
M_.param_names(35) = {'rho_z'};
M_.param_names_tex(35) = {'rho\_z'};
M_.param_names_long(35) = {'rho_z'};
M_.param_names(36) = {'sigma_z'};
M_.param_names_tex(36) = {'sigma\_z'};
M_.param_names_long(36) = {'sigma_z'};
M_.param_names(37) = {'rho_g'};
M_.param_names_tex(37) = {'rho\_g'};
M_.param_names_long(37) = {'rho_g'};
M_.param_names(38) = {'sigma_g'};
M_.param_names_tex(38) = {'sigma\_g'};
M_.param_names_long(38) = {'sigma_g'};
M_.param_names(39) = {'c_ss'};
M_.param_names_tex(39) = {'c\_ss'};
M_.param_names_long(39) = {'c_ss'};
M_.param_names(40) = {'lambda_ss'};
M_.param_names_tex(40) = {'lambda\_ss'};
M_.param_names_long(40) = {'lambda_ss'};
M_.param_names(41) = {'h_ss'};
M_.param_names_tex(41) = {'h\_ss'};
M_.param_names_long(41) = {'h_ss'};
M_.param_names(42) = {'l_ss'};
M_.param_names_tex(42) = {'l\_ss'};
M_.param_names_long(42) = {'l_ss'};
M_.param_names(43) = {'w_ss'};
M_.param_names_tex(43) = {'w\_ss'};
M_.param_names_long(43) = {'w_ss'};
M_.param_names(44) = {'R_ss'};
M_.param_names_tex(44) = {'R\_ss'};
M_.param_names_long(44) = {'R_ss'};
M_.param_names(45) = {'pi_ss'};
M_.param_names_tex(45) = {'pi\_ss'};
M_.param_names_long(45) = {'pi_ss'};
M_.param_names(46) = {'q_ss'};
M_.param_names_tex(46) = {'q\_ss'};
M_.param_names_long(46) = {'q_ss'};
M_.param_names(47) = {'r_k_ss'};
M_.param_names_tex(47) = {'r\_k\_ss'};
M_.param_names_long(47) = {'r_k_ss'};
M_.param_names(48) = {'Nfy_ss'};
M_.param_names_tex(48) = {'Nfy\_ss'};
M_.param_names_long(48) = {'Nfy_ss'};
M_.param_names(49) = {'ukl_ss'};
M_.param_names_tex(49) = {'ukl\_ss'};
M_.param_names_long(49) = {'ukl_ss'};
M_.param_names(50) = {'gy_ss'};
M_.param_names_tex(50) = {'gy\_ss'};
M_.param_names_long(50) = {'gy_ss'};
M_.param_names(51) = {'ky_ss'};
M_.param_names_tex(51) = {'ky\_ss'};
M_.param_names_long(51) = {'ky_ss'};
M_.param_names(52) = {'xk_ss'};
M_.param_names_tex(52) = {'xk\_ss'};
M_.param_names_long(52) = {'xk_ss'};
M_.param_names(53) = {'k_ss'};
M_.param_names_tex(53) = {'k\_ss'};
M_.param_names_long(53) = {'k_ss'};
M_.param_names(54) = {'x_ss'};
M_.param_names_tex(54) = {'x\_ss'};
M_.param_names_long(54) = {'x_ss'};
M_.param_names(55) = {'y_ss'};
M_.param_names_tex(55) = {'y\_ss'};
M_.param_names_long(55) = {'y_ss'};
M_.param_names(56) = {'N_ss'};
M_.param_names_tex(56) = {'N\_ss'};
M_.param_names_long(56) = {'N_ss'};
M_.param_names(57) = {'u_ss'};
M_.param_names_tex(57) = {'u\_ss'};
M_.param_names_long(57) = {'u_ss'};
M_.param_names(58) = {'N_star_ss'};
M_.param_names_tex(58) = {'N\_star\_ss'};
M_.param_names_long(58) = {'N_star_ss'};
M_.param_names(59) = {'Omega_ss'};
M_.param_names_tex(59) = {'Omega\_ss'};
M_.param_names_long(59) = {'Omega_ss'};
M_.param_names(60) = {'theta_p_ss'};
M_.param_names_tex(60) = {'theta\_p\_ss'};
M_.param_names_long(60) = {'theta_p_ss'};
M_.param_names(61) = {'pi_w_ss'};
M_.param_names_tex(61) = {'pi\_w\_ss'};
M_.param_names_long(61) = {'pi_w_ss'};
M_.param_names(62) = {'theta_w_ss'};
M_.param_names_tex(62) = {'theta\_w\_ss'};
M_.param_names_long(62) = {'theta_w_ss'};
M_.param_names(63) = {'v_ss'};
M_.param_names_tex(63) = {'v\_ss'};
M_.param_names_long(63) = {'v_ss'};
M_.param_names(64) = {'xi_ss'};
M_.param_names_tex(64) = {'xi\_ss'};
M_.param_names_long(64) = {'xi_ss'};
M_.param_names(65) = {'phi_ss'};
M_.param_names_tex(65) = {'phi\_ss'};
M_.param_names_long(65) = {'phi_ss'};
M_.param_names(66) = {'psi_ss'};
M_.param_names_tex(66) = {'psi\_ss'};
M_.param_names_long(66) = {'psi_ss'};
M_.param_names(67) = {'z_ss'};
M_.param_names_tex(67) = {'z\_ss'};
M_.param_names_long(67) = {'z_ss'};
M_.param_names(68) = {'g_ss'};
M_.param_names_tex(68) = {'g\_ss'};
M_.param_names_long(68) = {'g_ss'};
M_.param_names(69) = {'sig_p'};
M_.param_names_tex(69) = {'sig\_p'};
M_.param_names_long(69) = {'sig_p'};
M_.param_names(70) = {'sig_w'};
M_.param_names_tex(70) = {'sig\_w'};
M_.param_names_long(70) = {'sig_w'};
M_.param_names(71) = {'sig_n'};
M_.param_names_tex(71) = {'sig\_n'};
M_.param_names_long(71) = {'sig_n'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 10;
M_.endo_nbr = 26;
M_.param_nbr = 71;
M_.orig_endo_nbr = 26;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(10, 10);
M_.Correlation_matrix = eye(10, 10);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = true;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 26;
M_.eq_nbr = 26;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 1 15 41;
 2 16 0;
 0 17 0;
 3 18 0;
 4 19 0;
 5 20 42;
 0 21 43;
 0 22 44;
 6 23 0;
 7 24 45;
 0 25 0;
 8 26 0;
 0 27 0;
 0 28 0;
 0 29 0;
 0 30 46;
 0 31 0;
 9 32 0;
 10 33 0;
 11 34 0;
 12 35 0;
 13 36 0;
 14 37 0;
 0 38 0;
 0 39 0;
 0 40 0;]';
M_.nstatic = 9;
M_.nfwrd   = 3;
M_.npred   = 11;
M_.nboth   = 3;
M_.nsfwrd   = 6;
M_.nspred   = 14;
M_.ndynamic   = 17;
M_.dynamic_tmp_nbr = [8; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'c' ;
  2 , 'name' , 'x' ;
  3 , 'name' , '3' ;
  4 , 'name' , 'h' ;
  5 , 'name' , 'q' ;
  6 , 'name' , 'k' ;
  7 , 'name' , 'y' ;
  8 , 'name' , 'N' ;
  9 , 'name' , '9' ;
  10 , 'name' , 'u' ;
  11 , 'name' , 'pi' ;
  12 , 'name' , 'theta_p' ;
  13 , 'name' , '13' ;
  14 , 'name' , '14' ;
  15 , 'name' , 'theta_w' ;
  16 , 'name' , '16' ;
  17 , 'name' , 'R' ;
  18 , 'name' , 'v' ;
  19 , 'name' , 'xi' ;
  20 , 'name' , 'phi' ;
  21 , 'name' , 'psi' ;
  22 , 'name' , 'z' ;
  23 , 'name' , 'g' ;
  24 , 'name' , 'eps_p' ;
  25 , 'name' , 'eps_w' ;
  26 , 'name' , 'eps_n' ;
};
M_.mapping.c.eqidx = [1 4 16 ];
M_.mapping.h.eqidx = [3 4 ];
M_.mapping.l.eqidx = [3 7 9 14 ];
M_.mapping.w.eqidx = [3 13 14 ];
M_.mapping.R.eqidx = [1 5 17 ];
M_.mapping.pi.eqidx = [1 5 11 13 14 17 ];
M_.mapping.q.eqidx = [2 5 ];
M_.mapping.r_k.eqidx = [5 9 10 ];
M_.mapping.k.eqidx = [6 7 9 14 ];
M_.mapping.x.eqidx = [2 6 16 ];
M_.mapping.y.eqidx = [7 8 16 17 ];
M_.mapping.N.eqidx = [7 8 9 12 14 15 ];
M_.mapping.u.eqidx = [7 9 10 14 ];
M_.mapping.Omega.eqidx = [9 11 14 ];
M_.mapping.theta_p.eqidx = [8 11 12 ];
M_.mapping.pi_w.eqidx = [13 14 ];
M_.mapping.theta_w.eqidx = [8 14 15 ];
M_.mapping.v.eqidx = [1 18 ];
M_.mapping.xi.eqidx = [3 19 ];
M_.mapping.phi.eqidx = [1 5 20 ];
M_.mapping.psi.eqidx = [2 6 21 ];
M_.mapping.z.eqidx = [7 9 14 22 ];
M_.mapping.g.eqidx = [16 23 ];
M_.mapping.eps_p.eqidx = [11 24 ];
M_.mapping.eps_w.eqidx = [14 25 ];
M_.mapping.eps_n.eqidx = [8 26 ];
M_.mapping.e_v.eqidx = [18 ];
M_.mapping.e_xi.eqidx = [19 ];
M_.mapping.e_phi.eqidx = [20 ];
M_.mapping.e_psi.eqidx = [21 ];
M_.mapping.e_z.eqidx = [22 ];
M_.mapping.e_R.eqidx = [17 ];
M_.mapping.e_g.eqidx = [23 ];
M_.mapping.e_p.eqidx = [24 ];
M_.mapping.e_w.eqidx = [25 ];
M_.mapping.e_n.eqidx = [26 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 2 4 5 6 9 10 12 18 19 20 21 22 23 ];
M_.exo_names_orig_ord = [1:10];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(26, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(10, 1);
M_.params = NaN(71, 1);
M_.endo_trends = struct('deflator', cell(26, 1), 'log_deflator', cell(26, 1), 'growth_factor', cell(26, 1), 'log_growth_factor', cell(26, 1));
M_.NNZDerivatives = [114; 0; -1; ];
M_.static_tmp_nbr = [9; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 1.005;
gamma = M_.params(1);
M_.params(4) = 0.995;
beta = M_.params(4);
M_.params(5) = 0.02;
delta = M_.params(5);
M_.params(8) = 0.3;
alpha = M_.params(8);
M_.params(13) = 0.45;
s_m = M_.params(13);
M_.params(50) = 0.2;
gy_ss = M_.params(50);
M_.params(42) = 3;
l_ss = M_.params(42);
M_.params(56) = 9;
N_ss = M_.params(56);
M_.params(2) = 0.9736;
zeta = M_.params(2);
M_.params(16) = 9.9582;
eta_y = M_.params(16);
M_.params(15) = 1.9790;
chi_y = M_.params(15);
M_.params(18) = 11.1791;
chi_l = M_.params(18);
M_.params(3) = 2.5078;
vartheta = M_.params(3);
M_.params(19) = 1/M_.params(3);
eta_l = M_.params(19);
M_.params(7) = 0.003;
rho = M_.params(7);
M_.params(6) = 6.5137;
kappa_x = M_.params(6);
M_.params(12) = 0.3053;
omega = M_.params(12);
M_.params(14) = 236.3397;
kappa_p = M_.params(14);
M_.params(17) = 168.8219;
kappa_w = M_.params(17);
M_.params(10) = 0.9737;
rho_N = M_.params(10);
M_.params(21) = 0.9015;
gamma_p = M_.params(21);
M_.params(22) = 0.6419;
gamma_w = M_.params(22);
M_.params(63) = 1;
v_ss = M_.params(63);
M_.params(64) = 1;
xi_ss = M_.params(64);
M_.params(65) = 1;
phi_ss = M_.params(65);
M_.params(66) = 1;
psi_ss = M_.params(66);
M_.params(67) = 1;
z_ss = M_.params(67);
M_.params(57) = 1;
u_ss = M_.params(57);
M_.params(58) = M_.params(56);
N_star_ss = M_.params(58);
M_.params(45) = 1.005;
pi_ss = M_.params(45);
M_.params(61) = M_.params(45);
pi_w_ss = M_.params(61);
M_.params(46) = 1;
q_ss = M_.params(46);
M_.params(44) = M_.params(45)/M_.params(4);
R_ss = M_.params(44);
M_.params(60) = (1-M_.params(13))*(M_.params(15)-(M_.params(15)-M_.params(16))/M_.params(56))/((1-M_.params(13))*(M_.params(15)-(M_.params(15)-M_.params(16))/M_.params(56))-1);
theta_p_ss = M_.params(60);
M_.params(62) = (M_.params(18)-(M_.params(18)-M_.params(19))/M_.params(56))/(1+M_.params(18)-(M_.params(18)-M_.params(19))/M_.params(56));
theta_w_ss = M_.params(62);
M_.params(48) = M_.params(60)/(M_.params(8)+(1-M_.params(8))*M_.params(62))-1;
Nfy_ss = M_.params(48);
M_.params(59) = 1/M_.params(60);
Omega_ss = M_.params(59);
M_.params(52) = M_.params(1)-1+M_.params(5);
xk_ss = M_.params(52);
M_.params(47) = M_.params(46)*(M_.params(5)+M_.params(1)/M_.params(4)-1);
r_k_ss = M_.params(47);
M_.params(49) = (M_.params(47)/M_.params(59)/M_.params(8)/M_.params(67))^(1/(M_.params(8)-1));
ukl_ss = M_.params(49);
M_.params(43) = M_.params(67)*(1-M_.params(8))*M_.params(62)*M_.params(59)*M_.params(49)^M_.params(8);
w_ss = M_.params(43);
M_.params(51) = M_.params(8)*(1+M_.params(48))/(M_.params(60)*(M_.params(5)+1/M_.params(4)-1));
ky_ss = M_.params(51);
M_.params(53) = M_.params(49)/M_.params(57)*M_.params(42);
k_ss = M_.params(53);
M_.params(54) = M_.params(52)*M_.params(53);
x_ss = M_.params(54);
M_.params(55) = M_.params(53)/M_.params(51);
y_ss = M_.params(55);
M_.params(68) = M_.params(50)*M_.params(55);
g_ss = M_.params(68);
M_.params(39) = M_.params(55)-M_.params(54)-M_.params(68);
c_ss = M_.params(39);
M_.params(9) = M_.params(55)*M_.params(48)/M_.params(56);
f = M_.params(9);
M_.params(41) = (1-M_.params(2))*M_.params(39);
h_ss = M_.params(41);
M_.params(40) = M_.params(63)/((1-M_.params(2))*M_.params(39));
lambda_ss = M_.params(40);
M_.params(11) = M_.params(8)*M_.params(59)*(M_.params(55)/M_.params(56)+M_.params(9))/(M_.params(53)*M_.params(57)^(1+M_.params(12))/M_.params(56));
kappa_u = M_.params(11);
M_.params(20) = M_.params(43)/(M_.params(64)*M_.params(41)*M_.params(42)^M_.params(3));
gamma_l = M_.params(20);
M_.params(23) = 0.8715;
rho_R = M_.params(23);
M_.params(24) = 1.7946;
alpha_pi = M_.params(24);
M_.params(25) = 0.2065;
alpha_y = M_.params(25);
M_.params(26) = 0.0009;
sigma_R = M_.params(26);
M_.params(27) = 0.1192;
rho_v = M_.params(27);
M_.params(28) = 0.0019;
sigma_v = M_.params(28);
M_.params(29) = 0.8472;
rho_xi = M_.params(29);
M_.params(30) = 0.0188;
sigma_xi = M_.params(30);
M_.params(31) = 0.7450;
rho_phi = M_.params(31);
M_.params(32) = 0.0198;
sigma_phi = M_.params(32);
M_.params(33) = 0.4980;
rho_psi = M_.params(33);
M_.params(34) = 0.0023;
sigma_psi = M_.params(34);
M_.params(35) = 0.8266;
rho_z = M_.params(35);
M_.params(36) = 0.0044;
sigma_z = M_.params(36);
M_.params(37) = 0.8514;
rho_g = M_.params(37);
M_.params(38) = 0.0151;
sigma_g = M_.params(38);
M_.params(69) = .0022;
sig_p = M_.params(69);
M_.params(70) = .0101;
sig_w = M_.params(70);
M_.params(71) = .0038;
sig_n = M_.params(71);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(1) = 0;
oo_.steady_state(2) = 0;
oo_.steady_state(3) = 0;
oo_.steady_state(4) = 0;
oo_.steady_state(5) = 0;
oo_.steady_state(6) = 0;
oo_.steady_state(7) = 0;
oo_.steady_state(8) = 0;
oo_.steady_state(9) = 0;
oo_.steady_state(10) = 0;
oo_.steady_state(11) = 0;
oo_.steady_state(12) = 0;
oo_.steady_state(13) = 0;
oo_.steady_state(14) = 0;
oo_.steady_state(15) = 0;
oo_.steady_state(16) = 0;
oo_.steady_state(17) = 0;
oo_.steady_state(18) = 0;
oo_.steady_state(19) = 0;
oo_.steady_state(20) = 0;
oo_.steady_state(21) = 0;
oo_.steady_state(22) = 0;
oo_.steady_state(23) = 0;
oo_.steady_state(24) = 0;
oo_.steady_state(26) = 0;
oo_.steady_state(25) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
resid(1);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
model_diagnostics(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 0;
M_.Sigma_e(2, 2) = 0;
M_.Sigma_e(3, 3) = 0;
M_.Sigma_e(4, 4) = 0;
M_.Sigma_e(5, 5) = 1;
M_.Sigma_e(6, 6) = 1;
M_.Sigma_e(7, 7) = 0;
M_.Sigma_e(8, 8) = 0;
M_.Sigma_e(9, 9) = 1;
M_.Sigma_e(10, 10) = 0;
options_.irf = 20;
options_.order = 1;
var_list_ = {'y';'c';'l';'w';'pi';'theta_p';'theta_w';'N';'r_k';'R'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'kang_linear_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'kang_linear_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'kang_linear_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'kang_linear_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'kang_linear_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'kang_linear_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'kang_linear_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
