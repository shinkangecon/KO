function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = kang_linear.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(26, 56);
g1(1,1)=(-(T(1)/(1+T(1))));
g1(1,15)=1;
g1(1,41)=(-(1/(1+T(1))));
g1(1,19)=T(2);
g1(1,42)=(-T(2));
g1(1,32)=(-1);
g1(1,34)=T(2);
g1(2,21)=(-T(3));
g1(2,7)=(-(1/(1+params(4))));
g1(2,24)=1;
g1(2,45)=(-(params(4)/(1+params(4))));
g1(2,35)=(-T(3));
g1(3,16)=1;
g1(3,17)=params(3);
g1(3,18)=(-1);
g1(3,33)=1;
g1(4,1)=(-(params(7)/(1-T(1))*(-T(1))));
g1(4,15)=(-(params(7)/(1-T(1))));
g1(4,2)=(-((1-params(7))/params(1)));
g1(4,16)=1;
g1(5,19)=1;
g1(5,42)=(-1);
g1(5,21)=1;
g1(5,43)=(-(params(4)*(1-params(5))/params(1)));
g1(5,44)=(-(1-params(4)*(1-params(5))/params(1)));
g1(5,34)=1;
g1(6,6)=(-((1-params(5))/params(1)));
g1(6,23)=1;
g1(6,24)=(-(1-(1-params(5))/params(1)));
g1(6,35)=(-(1-(1-params(5))/params(1)));
g1(7,17)=(-((1+params(48))*(1-params(8))));
g1(7,6)=(-((1+params(48))*params(8)));
g1(7,25)=1;
g1(7,8)=(-((1+params(48))*(-params(8))));
g1(7,26)=(-((1+params(48))*params(8)-params(48)));
g1(7,27)=(-((1+params(48))*params(8)));
g1(7,36)=(-(1+params(48)));
g1(8,25)=(-(1-params(10)));
g1(8,8)=(-params(10));
g1(8,26)=1;
g1(8,29)=(-((1-params(10))*T(4)));
g1(8,31)=(-((1-params(10))*T(4)*(-((1-params(8))*params(62)/(params(8)+(1-params(8))*params(62))))));
g1(8,40)=(-1);
g1(9,17)=(-(params(8)-1));
g1(9,22)=(-1);
g1(9,6)=params(8)-1;
g1(9,8)=(-(params(8)-1));
g1(9,26)=params(8)-1;
g1(9,27)=params(8)-1;
g1(9,28)=1;
g1(9,36)=1;
g1(10,22)=(-(1/params(12)));
g1(10,27)=1;
g1(11,5)=(-(params(21)/(1+params(4)*params(21))));
g1(11,20)=1;
g1(11,42)=(-(params(4)/(1+params(4)*params(21))));
g1(11,28)=(-T(5));
g1(11,29)=(-T(5));
g1(11,38)=(-1);
g1(12,26)=(-T(6));
g1(12,29)=1;
g1(13,3)=1;
g1(13,18)=(-1);
g1(13,20)=(-1);
g1(13,30)=1;
g1(14,17)=(-(T(7)*(-params(8))));
g1(14,18)=T(7);
g1(14,5)=(-params(22));
g1(14,20)=(-(params(4)*(-params(22))));
g1(14,6)=(-(params(8)*T(7)));
g1(14,8)=(-(T(7)*(-params(8))));
g1(14,26)=(-(params(8)*T(7)));
g1(14,27)=(-(params(8)*T(7)));
g1(14,28)=(-T(7));
g1(14,30)=1;
g1(14,46)=(-params(4));
g1(14,31)=(-T(7));
g1(14,36)=(-T(7));
g1(14,39)=(-1);
g1(15,26)=(-T(8));
g1(15,31)=1;
g1(16,15)=params(39);
g1(16,24)=params(54);
g1(16,25)=(-params(55));
g1(16,37)=params(68);
g1(17,4)=(-params(23));
g1(17,19)=1;
g1(17,20)=(-((1-params(23))*params(24)));
g1(17,25)=(-((1-params(23))*(-params(25))));
g1(17,52)=(-params(26));
g1(18,9)=(-params(27));
g1(18,32)=1;
g1(18,47)=(-params(28));
g1(19,10)=(-params(29));
g1(19,33)=1;
g1(19,48)=(-params(30));
g1(20,11)=(-params(31));
g1(20,34)=1;
g1(20,49)=(-params(32));
g1(21,12)=(-params(33));
g1(21,35)=1;
g1(21,50)=(-params(34));
g1(22,13)=(-params(35));
g1(22,36)=1;
g1(22,51)=(-params(36));
g1(23,14)=(-params(37));
g1(23,37)=1;
g1(23,53)=(-params(38));
g1(24,38)=1;
g1(24,54)=(-params(69));
g1(25,39)=1;
g1(25,55)=(-params(70));
g1(26,40)=1;
g1(26,56)=(-params(71));

end
