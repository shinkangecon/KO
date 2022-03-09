function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = kang_linear.static_g1_tt(T, y, x, params);
end
g1 = zeros(26, 26);
g1(1,1)=1-(T(1)/(1+T(1))+1/(1+T(1)));
g1(1,5)=T(2);
g1(1,6)=(-T(2));
g1(1,18)=(-1);
g1(1,20)=T(2);
g1(2,7)=(-T(3));
g1(2,10)=1-(1/(1+params(4))+params(4)/(1+params(4)));
g1(2,21)=(-T(3));
g1(3,2)=1;
g1(3,3)=params(3);
g1(3,4)=(-1);
g1(3,19)=1;
g1(4,1)=(-params(7));
g1(4,2)=1-(1-params(7))/params(1);
g1(5,5)=1;
g1(5,6)=(-1);
g1(5,7)=T(4);
g1(5,8)=(-T(4));
g1(5,20)=1;
g1(6,9)=1-(1-params(5))/params(1);
g1(6,10)=(-(1-(1-params(5))/params(1)));
g1(6,21)=(-(1-(1-params(5))/params(1)));
g1(7,3)=(-((1+params(48))*(1-params(8))));
g1(7,9)=(-((1+params(48))*params(8)));
g1(7,11)=1;
g1(7,12)=params(48);
g1(7,13)=(-((1+params(48))*params(8)));
g1(7,22)=(-(1+params(48)));
g1(8,11)=(-(1-params(10)));
g1(8,12)=1-params(10);
g1(8,15)=(-((1-params(10))*T(5)));
g1(8,17)=(-((1-params(10))*T(5)*(-((1-params(8))*params(62)/(params(8)+(1-params(8))*params(62))))));
g1(8,26)=(-1);
g1(9,3)=(-(params(8)-1));
g1(9,8)=(-1);
g1(9,9)=params(8)-1;
g1(9,13)=params(8)-1;
g1(9,14)=1;
g1(9,22)=1;
g1(10,8)=(-(1/params(12)));
g1(10,13)=1;
g1(11,6)=1-(params(21)/(1+params(4)*params(21))+params(4)/(1+params(4)*params(21)));
g1(11,14)=(-T(6));
g1(11,15)=(-T(6));
g1(11,24)=(-1);
g1(12,12)=(-T(7));
g1(12,15)=1;
g1(13,6)=(-1);
g1(13,16)=1;
g1(14,3)=(-(T(8)*(-params(8))));
g1(14,4)=T(8);
g1(14,6)=(-params(22))-params(4)*(-params(22));
g1(14,9)=(-(params(8)*T(8)));
g1(14,13)=(-(params(8)*T(8)));
g1(14,14)=(-T(8));
g1(14,16)=1-params(4);
g1(14,17)=(-T(8));
g1(14,22)=(-T(8));
g1(14,25)=(-1);
g1(15,12)=(-T(9));
g1(15,17)=1;
g1(16,1)=params(39);
g1(16,10)=params(54);
g1(16,11)=(-params(55));
g1(16,23)=params(68);
g1(17,5)=1-params(23);
g1(17,6)=(-((1-params(23))*params(24)));
g1(17,11)=(-((1-params(23))*(-params(25))));
g1(18,18)=1-params(27);
g1(19,19)=1-params(29);
g1(20,20)=1-params(31);
g1(21,21)=1-params(33);
g1(22,22)=1-params(35);
g1(23,23)=1-params(37);
g1(24,24)=1;
g1(25,25)=1;
g1(26,26)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
