clear
close all


%% Initializing Flight Parameters

%%%
% Atmospheric & Flight Conditions
%%%

g           = 32.17;                % ft / s^2

altitude    = 0;                    % ft
M           = 0.158;                % n/a

rho         = 2.3769 * 10 ^ (-3);   % slugs / ft^3
a           = 1116.45;              % ft / s

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

C_L             =  0.41; 
C_D             =  0.054;

C_L_alpha       =  4.44;
C_D_alpha       =  0.33;
C_m_alpha       = -0.683;

C_L_alpha_dot   =  0;
C_m_alpha_dot   = -4.36;

C_L_q           =  3.8;
C_m_q           = -9.96;

C_L_M           =  0.0;
C_D_M           =  0.0;
C_m_M           =  0.0;

C_L_delta_e     =  0.355;
C_m_delta_e     = -0.923;


%%%
% Aircraft Lateral Coefficients (all derivatives are per radian)
%%%

C_y_beta        = -0.564;
C_l_beta        = -0.074;
C_n_beta        =  0.071;

C_l_p           = -0.410;
C_n_p           = -0.0575;

C_l_r           =  0.107;
C_n_r           = -0.125;

C_l_delta_a     = -0.134;
C_n_delta_a     =  0.0035;

C_y_delta_r     =  0.157;
C_l_delta_r     =  0.107;
C_n_delta_r     = -0.072;

%%%
% Aircraft Center of Gravity & Mass Characteristics
%%%

W       = 2750;     % lb
m       = W / g;    % slugs
I_x     = 1048;     % slugs*ft^2
I_y     = 3000;     % slugs*ft^2
I_z     = 3530;     % slugs*ft^2
I_xz    =    0;     % slugs*ft^2

%%%
% Aircraft Reference Geometry
%%%

S       = 184;      % ft^2
b       =  33.4;    % ft
c_bar   =   5.7;    % ft

