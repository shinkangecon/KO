function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 8);

T(1) = params(2)/params(1);
T(2) = (1-T(1))/(1+T(1));
T(3) = 1/(1+params(4))/params(6);
T(4) = params(60)/(params(60)-(params(8)+(1-params(8))*params(62)));
T(5) = ((1-params(13))*(params(15)-(params(15)-params(16))/params(56))-1)/(1+params(4)*params(21))/params(14);
T(6) = (params(15)-params(16))*(-(1-params(13)))/params(56)*(params(60)-1)^2/params(60);
T(7) = (1+params(18)-(params(18)-params(19))/params(56))/params(17);
T(8) = (params(18)-params(19))/params(56)*(1-params(62))^2/params(62);

end
