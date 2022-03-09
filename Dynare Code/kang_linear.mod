%% Business Cycle Implications of Firm Market Power in Labor and Product Markets
%% Alpanda and Zubairy (2021)
%% Kang and Oh JME GAZUAHHH!!

var c h l w R pi q r_k k x y N u Omega theta_p pi_w theta_w % #19
    v xi phi psi z g eps_p eps_w eps_n% #9
    ; 

%predetermined_variables N;




varexo e_v e_xi e_phi e_psi e_z e_R e_g e_p e_w e_n % #7
        ;

parameters gamma zeta vartheta beta delta kappa_x rho alpha f rho_N kappa_u omega s_m kappa_p chi_y eta_y kappa_w chi_l eta_l gamma_l gamma_p gamma_w
        rho_R alpha_pi alpha_y sigma_R rho_v sigma_v rho_xi sigma_xi rho_phi sigma_phi rho_psi sigma_psi rho_z sigma_z rho_g sigma_g
        c_ss lambda_ss h_ss l_ss w_ss R_ss pi_ss q_ss r_k_ss Nfy_ss ukl_ss gy_ss ky_ss xk_ss k_ss x_ss y_ss N_ss u_ss N_star_ss Omega_ss theta_p_ss pi_w_ss theta_w_ss
        v_ss xi_ss phi_ss psi_ss z_ss g_ss sig_p sig_w sig_n;


gamma = 1.005;
beta = 0.995;
delta = 0.02;
alpha = 0.3;
s_m = 0.45;
gy_ss = 0.2;
l_ss = 3; % s.t. gamma_l = 1
N_ss = 9;
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
xk_ss = gamma - 1 + delta;
r_k_ss = q_ss*(gamma/beta - 1 + delta);
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
gamma_l = w_ss/(xi_ss*h_ss*l_ss^vartheta);

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
sig_p=.0022;
sig_w=.0101;
sig_n=.0038;
 
model(linear);
% 1.
    c = zeta/gamma/(1 + zeta/gamma)*c(-1) + 1/(1 + zeta/gamma)*c(+1) - (1 - zeta/gamma)/(1 + zeta/gamma)*(R - pi(+1) + phi) + v;

% 2.
    x = 1/(1 + beta)*x(-1) + beta/(1 + beta)*x(+1) + 1/(1 + beta)/kappa_x*(q + psi);

% 3.
    xi + h + vartheta*l = w;

% 4. 
    h = rho/(1 - zeta/gamma)*(c - zeta/gamma*c(-1)) + (1 - rho)/gamma*h(-1);

% 5.
    q = (1 - delta)*beta/gamma*q(+1) + (1 - (1 - delta)*beta/gamma)*r_k(+1) - (R - pi(+1) + phi);

% 6. 1-delta/gamma->1+delta/gamma
    k = (1 - delta)/gamma*k(-1) + (1 - (1 - delta)/gamma)*(x + psi);

% 7.
    y = (1 + Nfy_ss)*(z + alpha*(u + N - N(-1) + k(-1)) + (1 - alpha)*l) - Nfy_ss*N;

% 8. eps_Nt?
    N = rho_N*N(-1) + (1 - rho_N)*(y + theta_p_ss/(theta_p_ss - ((1 - alpha)*theta_w_ss + alpha))*(theta_p - (1 - alpha)*theta_w_ss/((1 - alpha)*theta_w_ss + alpha)*theta_w))+eps_n;

% 9. 
    Omega + z + (alpha - 1)*(u + N - N(-1) + k(-1) - l) = r_k;

% 10.
    u = 1/omega*r_k;

% 11. eps_pt?
    pi = gamma_p/(1 + beta*gamma_p)*pi(-1) + beta/(1 + beta*gamma_p)*pi(+1) + ((1 - s_m)*(chi_y - (chi_y - eta_y)/N_ss) - 1)/(1 + beta*gamma_p)/kappa_p*(Omega + theta_p)+eps_p;

% 12. 
    theta_p = - (1 - s_m)*(chi_y - eta_y)/N_ss*(theta_p_ss - 1)^2/theta_p_ss*N;

% 13.
    pi_w - pi = w - w(-1);

% 14. eps_wt?
    pi_w - gamma_w*pi(-1) = beta*(pi_w(+1) - gamma_w*pi) + (chi_l - (chi_l - eta_l)/N_ss + 1)/kappa_w*(Omega + z + alpha*(u + N - N(-1) + k(-1) - l) - w + theta_w)+eps_w;

% 15.
    theta_w = (chi_l - eta_l)/N_ss*(1 - theta_w_ss)^2/theta_w_ss*N;

% 16.
    c_ss*c + x_ss*x + g_ss*g = y_ss*y;

% 17.
    R = rho_R*R(-1) + (1 - rho_R)*(alpha_pi*pi - alpha_y*y) + sigma_R*e_R;

% 20.
    (v) = rho_v*(v(-1)) + sigma_v*e_v;

% 21.
    (xi) = rho_xi*(xi(-1)) + sigma_xi*e_xi;

% 22.
    (phi) = rho_phi*(phi(-1)) + sigma_phi*e_phi;

% 23.
    (psi) = rho_psi*(psi(-1)) + sigma_psi*e_psi;

% 24.
    (z) = rho_z*(z(-1)) + sigma_z*e_z;

% 25.
    (g) = rho_g*(g(-1)) + sigma_g*e_g;
    
% 26.
    (eps_p)=sig_p*e_p;
    
% 27.
    (eps_w)=sig_w*e_w;

% 28.
    (eps_n)=sig_n*e_n;

end;

initval;
c = 0;
h = 0;
l = 0;
w = 0;
R = 0;
pi = 0;
q = 0;
r_k = 0;
k = 0;
x = 0;
y = 0;
N = 0;
u = 0;
Omega = 0;
theta_p = 0;
pi_w = 0;
theta_w = 0;
v = 0;
xi = 0;
phi = 0;
psi = 0;
z = 0;
g = 0;
eps_p=0;
eps_n=0;
eps_w=0;
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
var e_w=1;
var e_p=0;
var e_n=0;
end;

stoch_simul(order=1,irf=20) y c l w pi theta_p theta_w N r_k R;
