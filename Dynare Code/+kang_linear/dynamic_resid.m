function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = kang_linear.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(26, 1);
lhs = y(15);
rhs = T(1)/(1+T(1))*y(1)+1/(1+T(1))*y(41)-T(2)*(y(19)-y(42)+y(34))+y(32);
residual(1) = lhs - rhs;
lhs = y(24);
rhs = 1/(1+params(4))*y(7)+params(4)/(1+params(4))*y(45)+T(3)*(y(21)+y(35));
residual(2) = lhs - rhs;
lhs = y(33)+y(16)+params(3)*y(17);
rhs = y(18);
residual(3) = lhs - rhs;
lhs = y(16);
rhs = params(7)/(1-T(1))*(y(15)-T(1)*y(1))+(1-params(7))/params(1)*y(2);
residual(4) = lhs - rhs;
lhs = y(21);
rhs = params(4)*(1-params(5))/params(1)*y(43)+(1-params(4)*(1-params(5))/params(1))*y(44)-(y(19)-y(42)+y(34));
residual(5) = lhs - rhs;
lhs = y(23);
rhs = (1-params(5))/params(1)*y(6)+(1-(1-params(5))/params(1))*(y(24)+y(35));
residual(6) = lhs - rhs;
lhs = y(25);
rhs = (1+params(48))*(y(36)+params(8)*(y(6)+y(27)+y(26)-y(8))+y(17)*(1-params(8)))-params(48)*y(26);
residual(7) = lhs - rhs;
lhs = y(26);
rhs = y(8)*params(10)+(1-params(10))*(y(25)+T(4)*(y(29)-(1-params(8))*params(62)/(params(8)+(1-params(8))*params(62))*y(31)))+y(40);
residual(8) = lhs - rhs;
lhs = y(36)+y(28)+(params(8)-1)*(y(6)+y(27)+y(26)-y(8)-y(17));
rhs = y(22);
residual(9) = lhs - rhs;
lhs = y(27);
rhs = y(22)*1/params(12);
residual(10) = lhs - rhs;
lhs = y(20);
rhs = params(21)/(1+params(4)*params(21))*y(5)+y(42)*params(4)/(1+params(4)*params(21))+T(5)*(y(29)+y(28))+y(38);
residual(11) = lhs - rhs;
lhs = y(29);
rhs = y(26)*T(6);
residual(12) = lhs - rhs;
lhs = y(30)-y(20);
rhs = y(18)-y(3);
residual(13) = lhs - rhs;
lhs = y(30)-y(5)*params(22);
rhs = params(4)*(y(46)-y(20)*params(22))+T(7)*(y(31)+y(36)+y(28)+params(8)*(y(6)+y(27)+y(26)-y(8)-y(17))-y(18))+y(39);
residual(14) = lhs - rhs;
lhs = y(31);
rhs = y(26)*T(8);
residual(15) = lhs - rhs;
lhs = y(15)*params(39)+y(24)*params(54)+params(68)*y(37);
rhs = y(25)*params(55);
residual(16) = lhs - rhs;
lhs = y(19);
rhs = params(23)*y(4)+(1-params(23))*(y(20)*params(24)-y(25)*params(25))+params(26)*x(it_, 6);
residual(17) = lhs - rhs;
lhs = y(32);
rhs = params(27)*y(9)+params(28)*x(it_, 1);
residual(18) = lhs - rhs;
lhs = y(33);
rhs = params(29)*y(10)+params(30)*x(it_, 2);
residual(19) = lhs - rhs;
lhs = y(34);
rhs = params(31)*y(11)+params(32)*x(it_, 3);
residual(20) = lhs - rhs;
lhs = y(35);
rhs = params(33)*y(12)+params(34)*x(it_, 4);
residual(21) = lhs - rhs;
lhs = y(36);
rhs = params(35)*y(13)+params(36)*x(it_, 5);
residual(22) = lhs - rhs;
lhs = y(37);
rhs = params(37)*y(14)+params(38)*x(it_, 7);
residual(23) = lhs - rhs;
lhs = y(38);
rhs = params(69)*x(it_, 8);
residual(24) = lhs - rhs;
lhs = y(39);
rhs = params(70)*x(it_, 9);
residual(25) = lhs - rhs;
lhs = y(40);
rhs = params(71)*x(it_, 10);
residual(26) = lhs - rhs;

end
