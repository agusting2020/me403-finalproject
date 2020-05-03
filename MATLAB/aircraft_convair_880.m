clear
close all


%% Initializing Flight Parameters

%%%
% Atmospheric & Flight Conditions
%%%

g           = 32.17;                % ft / s^2

altitude    = 35000;                % ft
M           = 0.8;                  % n/a

rho         = 7.3820 * 10 ^ (-4);   % slugs / ft^3
a           = 973.14;               % ft / s

u_0         = M * a;                % ft / s
Q           = 0.5 * rho * u_0 ^ 2;  % psi


%%%
% Initial Conditions
%%%

x_0 = [
    0;
    0;
    1;
    5 * pi /180];

theta_0 = 0;

    
%% Initializing Aircraft Parameters

%%%
% Aircraft Longitudinal Coefficients (all derivatives are per radian)
%%%

C_L             =  0.347;
C_D             =  0.024;

C_L_alpha       =  4.8;
C_D_alpha       =  0.15;
C_m_alpha       = -0.65;

C_L_alpha_dot   =  2.7;
C_m_alpha_dot   = -4.5;

C_L_q           =  7.5;
C_m_q           = -4.5;

C_L_M           =  0.0;
C_D_M           =  0.0;
C_m_M           =  0.0;

C_L_delta_e     =  0.190;
C_m_delta_e     = -0.57;


%%%
% Aircraft Lateral Coefficients (all derivatives are per radian)
%%%

C_y_beta        = -0.812;
C_l_beta        = -0.177;
C_n_beta        =  0.129;

C_l_p           = -0.312;
C_n_p           = -0.011;

C_l_r           =  0.153;
C_n_r           = -0.165;

C_l_delta_a     = -0.050;
C_n_delta_a     =  0.008;

C_y_delta_r     =  0.184;
C_l_delta_r     =  0.019;
C_n_delta_r     = -0.076;

%%%
% Aircraft Center of Gravity & Mass Characteristics
%%%

W       = 126000;   % lb
m       = W / g;    % slugs
I_x     =  115000;  % slugs*ft^2
I_y     = 2450000;  % slugs*ft^2
I_z     = 4070000;  % slugs*ft^2
I_xz    =       0;  % slugs*ft^2

%%%
% Aircraft Reference Geometry
%%%

S       = 2000;     % ft^2
b       =  120;     % ft
c_bar   =   18.94;  % ft

