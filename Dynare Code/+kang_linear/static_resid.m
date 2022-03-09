function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = kang_linear.static_resid_tt(T, y, x, params);
end
residual = zeros(26, 1);
lhs = y(1);
rhs = y(1)*T(1)/(1+T(1))+y(1)*1/(1+T(1))-T(2)*(y(5)-y(6)+y(20))+y(18);
residual(1) = lhs - rhs;
lhs = y(10);
rhs = y(10)*1/(1+params(4))+y(10)*params(4)/(1+params(4))+T(3)*(y(7)+y(21));
residual(2) = lhs - rhs;
lhs = y(19)+y(2)+params(3)*y(3);
rhs = y(4);
residual(3) = lhs - rhs;
lhs = y(2);
rhs = params(7)/(1-T(1))*(y(1)-y(1)*T(1))+y(2)*(1-params(7))/params(1);
residual(4) = lhs - rhs;
lhs = y(7);
rhs = y(7)*params(4)*(1-params(5))/params(1)+T(4)*y(8)-(y(5)-y(6)+y(20));
residual(5) = lhs - rhs;
lhs = y(9);
rhs = y(9)*(1-params(5))/params(1)+(1-(1-params(5))/params(1))*(y(10)+y(21));
residual(6) = lhs - rhs;
lhs = y(11);
rhs = (1+params(48))*(y(22)+params(8)*(y(9)+y(13))+y(3)*(1-params(8)))-params(48)*y(12);
residual(7) = lhs - rhs;
lhs = y(12);
rhs = y(12)*params(10)+(1-params(10))*(y(11)+T(5)*(y(15)-(1-params(8))*params(62)/(params(8)+(1-params(8))*params(62))*y(17)))+y(26);
residual(8) = lhs - rhs;
lhs = y(22)+y(14)+(params(8)-1)*(y(9)+y(13)-y(3));
rhs = y(8);
residual(9) = lhs - rhs;
lhs = y(13);
rhs = y(8)*1/params(12);
residual(10) = lhs - rhs;
lhs = y(6);
rhs = y(6)*params(21)/(1+params(4)*params(21))+y(6)*params(4)/(1+params(4)*params(21))+T(6)*(y(15)+y(14))+y(24);
residual(11) = lhs - rhs;
lhs = y(15);
rhs = y(12)*T(7);
residual(12) = lhs - rhs;
residual(13) = y(16)-y(6);
lhs = y(16)-y(6)*params(22);
rhs = params(4)*(y(16)-y(6)*params(22))+T(8)*(y(17)+y(22)+y(14)+params(8)*(y(9)+y(13)-y(3))-y(4))+y(25);
residual(14) = lhs - rhs;
lhs = y(17);
rhs = y(12)*T(9);
residual(15) = lhs - rhs;
lhs = y(1)*params(39)+y(10)*params(54)+params(68)*y(23);
rhs = y(11)*params(55);
residual(16) = lhs - rhs;
lhs = y(5);
rhs = y(5)*params(23)+(1-params(23))*(y(6)*params(24)-y(11)*params(25))+params(26)*x(6);
residual(17) = lhs - rhs;
lhs = y(18);
rhs = y(18)*params(27)+params(28)*x(1);
residual(18) = lhs - rhs;
lhs = y(19);
rhs = y(19)*params(29)+params(30)*x(2);
residual(19) = lhs - rhs;
lhs = y(20);
rhs = y(20)*params(31)+params(32)*x(3);
residual(20) = lhs - rhs;
lhs = y(21);
rhs = y(21)*params(33)+params(34)*x(4);
residual(21) = lhs - rhs;
lhs = y(22);
rhs = y(22)*params(35)+params(36)*x(5);
residual(22) = lhs - rhs;
lhs = y(23);
rhs = y(23)*params(37)+params(38)*x(7);
residual(23) = lhs - rhs;
lhs = y(24);
rhs = params(69)*x(8);
residual(24) = lhs - rhs;
lhs = y(25);
rhs = params(70)*x(9);
residual(25) = lhs - rhs;
lhs = y(26);
rhs = params(71)*x(10);
residual(26) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
