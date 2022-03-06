%% Business Cycle Implications of Firm Market Power in Labor and Product Markets
%% Alpanda and Zubairy (2021)
%% Kang and Oh JME GAZUAHHH!!

var c lambda h l w R pi q r_k k x y N u N_star Omega theta_p pi_w theta_w % #19
    v xi phi psi z g% #5
    ; 

%predetermined_variables N;

varexo e_v e_xi e_phi e_psi e_z e_R e_g % #7
        ;

parameters zeta vartheta beta delta kappa_x rho alpha f rho_N kappa_u omega s_m kappa_p chi_y eta_y kappa_w chi_l eta_l gamma gamma_p gamma_w
        rho_R alpha_pi alpha_y sigma_R rho_v sigma_v rho_xi sigma_xi rho_phi sigma_phi rho_psi sigma_psi rho_z sigma_z rho_g sigma_g
        c_ss lambda_ss h_ss l_ss w_ss R_ss pi_ss q_ss r_k_ss ukl_ss gy_ss ky_ss xk_ss k_ss x_ss y_ss N_ss u_ss N_star_ss Omega_ss theta_p_ss pi_w_ss theta_w_ss
        v_ss xi_ss phi_ss psi_ss z_ss g_ss;


beta = 0.995;
delta = 0.02;
alpha = 0.3;
s_m = 0.45;
gy_ss = 0.2;
l_ss = 3; % s.t. gamma = 1
N_ss = 8;
zeta = 0.9736;
eta_y = 9.9582;
chi_y = 1.9790;
chi_l = 11.1791;
vartheta = 2.5078;
eta_l = 1/vartheta;
rho = 0.003;
kappa_x = 6.5137;
omega = 0.3053;
kappa_p = 236.3397;
kappa_w = 168.8219;
rho_N = 0.9737;
gamma_p = 0.9015;
gamma_w = 0.6419;
v_ss = 1;
xi_ss = 1;
phi_ss = 1;
psi_ss = 1;
z_ss = 1;
u_ss = 1;
N_star_ss = N_ss;
pi_ss = 1.005;
pi_w_ss = pi_ss;
q_ss = 1;
R_ss = pi_ss/beta;
theta_p_ss = (1 - s_m)*(chi_y - (chi_y - eta_y)/N_ss)/((1 - s_m)*(chi_y - (chi_y - eta_y)/N_ss) - 1);
theta_w_ss = (chi_l - (chi_l - eta_l)/N_ss)/(chi_l - (chi_l - eta_l)/N_ss + 1);
Nfy_ss = theta_p_ss/((1-alpha)*theta_w_ss + alpha) - 1;
Omega_ss = 1/theta_p_ss;
xk_ss = delta;
r_k_ss = q_ss*(1/beta - 1 + delta);
ukl_ss = (r_k_ss/Omega_ss/alpha/z_ss)^(1/(alpha-1));  
w_ss = theta_w_ss*Omega_ss*(1-alpha)*z_ss*ukl_ss^alpha;
ky_ss = alpha*(1 + Nfy_ss)/(theta_p_ss*(1/beta -1 + delta));
k_ss = ukl_ss/u_ss*l_ss;
x_ss = xk_ss*k_ss;
y_ss = k_ss/ky_ss;
g_ss = gy_ss*y_ss;
c_ss = y_ss - x_ss - g_ss;
f = Nfy_ss/N_ss*y_ss;
h_ss = ((1-zeta)*c_ss);
lambda_ss = v_ss/((1-zeta)*c_ss); 
kappa_u = Omega_ss*alpha*(y_ss/N_ss + f)/(u_ss^(omega+1)*k_ss/N_ss);
gamma = w_ss/(xi_ss*h_ss*l_ss^vartheta);

rho_R = 0.8715;
alpha_pi = 1.7946;
alpha_y = 0.2065;
sigma_R = 0.0009;
rho_v = 0.1192; 
sigma_v = 0.0019;
rho_xi = 0.8472;
sigma_xi = 0.0188;
rho_phi = 0.7450; 
sigma_phi = 0.0198;
rho_psi = 0.4980; 
sigma_psi = 0.0023;
rho_z = 0.8266;
sigma_z = 0.0044;
rho_g = 0.8514;
sigma_g = 0.0151;
 
model;
% 1.
    v/(c - zeta*c(-1)) = lambda;

% 2.
    gamma*xi*h*l^vartheta = w;

% 3.
    1 = beta*lambda(+1)/lambda*R*phi/pi(+1);

% 4.
    q = beta*lambda(+1)/lambda*((1 - delta)*q(+1) + r_k(+1));

% 5.
    k = (1 - delta)*k(-1) + (1 - kappa_x/2*(x/x(-1) - 1)^2)*psi*x;

% 6.
    1 = q*psi*(1 - kappa_x/2*(x/x(-1) - 1)^2 - kappa_x*(x/x(-1) - 1)*x/x(-1)) 
        + beta*lambda(+1)/lambda*q(+1)*psi(+1)*kappa_x*(x(+1)/x - 1)*(x(+1)/x)^2 ;

% 7.
    h = (c - zeta*c(-1))^rho*h(-1)^(1 - rho);

% 8.
    y/N = z*(u*k(-1)/N(-1))^alpha*(l/N)^(1 - alpha) - f;

% 9.
    N = N(-1)^rho_N*N_star^(1 - rho_N);

% 10.
    N_star = (theta_p/((1 - alpha)*theta_w + alpha) - 1)*y/f;

% 11.
    Omega*alpha*(y/N + f)/(k(-1)/N(-1)) = r_k + kappa_u/(1 + omega)*(u^(1 + omega) - 1);

% 12.
    r_k + kappa_u/(1 + omega)*(u^(1 + omega) - 1) = kappa_u*u^(1 + omega);

% 13.
    (pi/(pi(-1)^gamma_p*pi_ss^(1-gamma_p)) - 1)*pi/(pi(-1)^gamma_p*pi_ss^(1-gamma_p)) = beta*lambda(+1)/lambda*(pi(+1)/(pi^gamma_p*pi_ss^(1-gamma_p)) - 1)*pi(+1)/(pi(-1)^gamma_p*pi_ss^(1-gamma_p))*y(+1)/y*N/N(+1)
        - ((1 - s_m)*(chi_y - (chi_y - eta_y)/N) - 1)/kappa_p*(1 - theta_p*Omega);

% 14.
    theta_p = (1 - s_m)*(chi_y - (chi_y - eta_y)/N)/((1 - s_m)*(chi_y - (chi_y - eta_y)/N) - 1);

% 15.
    pi_w = w/w(-1)*pi;

% 16.
    (pi_w/(pi(-1)^gamma_w*pi_ss^(1-gamma_w)) - 1)*pi_w/(pi(-1)^gamma_w*pi_ss^(1-gamma_w)) = beta*lambda(+1)/lambda*(pi_w(+1)/(pi^gamma_w*pi_ss^(1-gamma_w)) - 1)*(pi_w(+1)/(pi^gamma_w*pi_ss^(1-gamma_w)))*(pi_w(+1)/pi_ss)*l(+1)/l*N/N(+1)
        - (chi_l - (chi_l - eta_l)/N + 1)/kappa_w*(1 - theta_w*Omega*(1-alpha)*z*(u*k(-1)*N/l/N(-1))^alpha/w);

% 17.
    theta_w = (chi_l - (chi_l - eta_l)/N)/(chi_l - (chi_l - eta_l)/N + 1);

% 18.
    c + x + g = y;

% 19.
    R/R_ss = (R(-1)/R)^rho_R*((pi/pi_ss)^alpha_pi*(y/y_ss)^alpha_y)^(1 - rho_R)*exp(sigma_R*e_R);

% 20.
    log(v) = rho_v*log(v(-1)) + sigma_v*e_v;

% 21.
    log(xi) = rho_xi*log(xi(-1)) + sigma_xi*e_xi;

% 22.
    log(phi) = rho_phi*log(phi(-1)) + sigma_phi*e_phi;

% 23.
    log(psi) = rho_psi*log(psi(-1)) + sigma_psi*e_psi;

% 24.
    log(z) = rho_z*log(z(-1)) + sigma_z*e_z;

% 25.
    log(g) = (1-rho_g)*log(g_ss) + rho_g*log(g(-1)) + sigma_g*e_g;

end;

initval;
c = c_ss;
lambda = lambda_ss; 
h = h_ss;
l = l_ss;
w = w_ss;
R = R_ss;
pi = pi_ss;
q = q_ss;
r_k = r_k_ss;
k = k_ss;
x = x_ss;
y = y_ss;
N = N_ss;
u = u_ss;
N_star = N_star_ss;
Omega = Omega_ss;
theta_p = theta_p_ss;
pi_w = pi_w_ss;
theta_w = theta_w_ss;
v = v_ss;
xi = xi_ss;
phi = phi_ss;
psi = psi_ss;
z = z_ss;
g = g_ss;
end;


resid(1);
steady;
check;

model_diagnostics;

shocks;
var e_v = 0;
var e_xi = 0;
var e_phi = 0;
var e_psi = 0;
var e_z = 1; 
var e_R = 1;
var e_g = 0;
end;

stoch_simul(order=1,irf=10) y c l w pi theta_p theta_w N r_k R;
