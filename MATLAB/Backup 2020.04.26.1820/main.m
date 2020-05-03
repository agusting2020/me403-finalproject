
%
%   Control Systems Analysis of the
%   Transport Aircraft: Convair 880
%   
%
%   Created by: Agustin Guzman
%


%% Missing Longitudinal Motion Coefficients

C_D_0           = C_D; %C_D = C_D_0;
C_L_0           = C_L; %C_L = C_L_0;

C_L_u           = 0;
C_D_u           = 0;
C_m_u           = 0;

C_z_alpha_dot   = 0;
C_Z_q           = 0;
C_Z_delta_e     = 0;
    

%% Missing Lateral Motion Coefficients

C_y_p       = 0;
C_y_r       = 0;
C_y_delta_a = 0;    


%% Initializing Derivatives

%%%
% Longitudinal Motion Derivatives
%%%

Z_u         = -(C_L_u + 2 * C_L_0) * Q * S / (m * u_0);
X_u         = -(C_D_u + 2 * C_D_0) * Q * S / (m * u_0);
M_u         = C_m_u * (Q * S * c_bar) / (u_0 * I_y);

Z_w         = -(C_L_alpha + C_D_0) * Q * S / (m * u_0);
X_w         = -(C_D_alpha - C_L_0) * Q * S / (m * u_0);
M_w         = C_m_alpha * (Q * S * c_bar) / (u_0 * I_y);

Z_w_dot     = C_z_alpha_dot * (c_bar / (2 * u_0)) * Q * S / (m * u_0);
M_w_dot     = C_m_alpha_dot * (c_bar / (2 * u_0)) * (Q * S * c_bar ) / (u_0 * I_y);

Z_alpha     = u_0 * Z_w;
M_alpha     = u_0 * M_w;

Z_alpha_dot = u_0 * Z_w_dot;
M_alpha_dot = u_0 * M_w_dot;

Z_q         = C_Z_q * (c_bar / (2 * u_0)) * Q * S / m;
M_q         = C_m_q * (c_bar / (2 * u_0)) * (Q * S * c_bar ) / I_y;

Z_delta_e   = C_Z_delta_e * Q * S / m;
M_delta_e   = C_m_delta_e * (Q * S * c_bar ) / I_y;


%%%
% Lateral Motion Derivatives
%%%

Y_beta      = Q * S * C_y_beta / m;
N_beta      = Q * S * b * C_n_beta / I_z;
L_beta      = Q * S * b * C_l_beta / I_x;

Y_p         = Q * S * b * C_y_p / (2 * m * u_0);
N_p         = Q * S * b ^ 2 * C_n_p / (2 * I_z * u_0);
L_p         = Q * S * b ^ 2 * C_l_p / (2 * I_x * u_0);

Y_r         = Q * S * b * C_y_r / (2 * m * u_0);
N_r         = Q * S * b ^ 2 * C_n_r / (2 * I_z * u_0);
L_r         = Q * S * b ^ 2 * C_l_r / (2 * I_x * u_0);

Y_delta_a   = Q * S * C_y_delta_a / m;
N_delta_a   = Q * S * b * C_n_delta_a / I_z;
L_delta_a   = Q * S * b * C_l_delta_a / I_x;

Y_delta_r   = Q * S * C_y_delta_r / m;
N_delta_r   = Q * S * b * C_n_delta_r / I_z;
L_delta_r   = Q * S * b * C_l_delta_r / I_x;

L_v         = L_beta;
N_v         = N_beta;
Y_v         = Y_beta;

L_star_v    = L_v / (1 - (I_xz ^ 2 / (I_x * I_z)));
N_star_v    = N_v / (1 - (I_xz ^ 2 / (I_x * I_z)));

L_star_r    = L_r / (1 - (I_xz ^ 2 / (I_x * I_z)));
N_star_r    = N_r / (1 - (I_xz ^ 2 / (I_x * I_z)));

L_star_p    = L_p / (1 - (I_xz ^ 2 / (I_x * I_z)));
N_star_p    = N_p / (1 - (I_xz ^ 2 / (I_x * I_z)));
    
    
%% Longitudinal Matrix Computation

%%%
% Longitudinal Motion Matrix A
%%%

Long_A_1_1 = X_u;
Long_A_1_2 = X_w;
Long_A_1_3 = 0;
Long_A_1_4 = -g;

Long_A_2_1 = Z_u;
Long_A_2_2 = Z_w;
Long_A_2_3 = u_0;
Long_A_2_4 = 0;

Long_A_3_1 = M_u + M_w_dot * Z_u;
Long_A_3_2 = M_w + M_w_dot * Z_w;
Long_A_3_3 = M_q + M_w_dot * u_0;
Long_A_3_4 = 0;

Long_A_4_1 = 0;
Long_A_4_2 = 0;
Long_A_4_3 = 1;
Long_A_4_4 = 0;

Long_A = [
    Long_A_1_1, Long_A_1_2, Long_A_1_3,	Long_A_1_4;
    Long_A_2_1, Long_A_2_2, Long_A_2_3, Long_A_2_4;
    Long_A_3_1, Long_A_3_2, Long_A_3_3, Long_A_3_4;
    Long_A_4_1, Long_A_4_2, Long_A_4_3, Long_A_4_4];

%%%
% Longitudinal Motion Matrix B
%%%

%     Long_B = [
%             X_delta,                                        X_delta_gamma;
%             Z_delta,                                        X_delta_gamma;
%             M_delta + M_w_dot * Z_delta,                    M_delta_gamma + M_w_dot * Z_delta_gamma;
% 
%             ];


%% Lateral Motion Matrix Computation

%%%
% Lateral Motion Matrix A
%%%

Lat_A_1_1 = Y_v;
Lat_A_1_2 = Y_p;
Lat_A_1_3 = -(u_0 - Y_r);
Lat_A_1_4 = g * cos(theta_0);

Lat_A_2_1 = L_star_v + (I_xz / I_x) * N_star_v;
Lat_A_2_2 = L_star_p + (I_xz / I_x) * N_star_p;
Lat_A_2_3 = L_star_r + (I_xz / I_x) * N_star_r;
Lat_A_2_4 = 0;

Lat_A_3_1 = N_star_v + (I_xz / I_z) * N_star_v;
Lat_A_3_2 = N_star_p + (I_xz / I_z) * N_star_p;
Lat_A_3_3 = N_star_r + (I_xz / I_z) * N_star_r;
Lat_A_3_4 = 0;

Lat_A_4_1 = 0;
Lat_A_4_2 = 1;
Lat_A_4_3 = 0;
Lat_A_4_4 = 0;

Lat_A = [
    Lat_A_1_1,   Lat_A_1_2,   Lat_A_1_3,   Lat_A_1_4;
    Lat_A_2_1,   Lat_A_2_2,   Lat_A_2_3,   Lat_A_2_4;
    Lat_A_3_1,   Lat_A_3_2,   Lat_A_3_3,   Lat_A_3_4;
    Lat_A_4_1,   Lat_A_4_2,   Lat_A_4_3,   Lat_A_4_4];


%%%
% Lateral Motion Matrix B
%%%

% Lat_B_1_1 = 0;
% Lat_B_1_2 = Y_delta_r;
% 
% Lat_B_2_1 = L_star_delta_a + (I_xz / I_x) * N_star_delta_a;
% Lat_B_2_2 = L_star_delta_r + (I_xz / I_x) * N_star_delta_r;
% 
% Lat_B_3_1 = N_star_delta_a + (I_xz / I_x) * L_star_delta_a;
% Lat_B_3_2 = N_star_delta_r + (I_xz / I_x) * L_star_delta_r;
% 
% Lat_B_4_1 = 0;
% Lat_B_4_2 = 0;
% 
% Lat_B = [
%         Lat_B_1_1,   Lat_B_1_2;
%         Lat_B_2_1,   Lat_B_2_2;
%         Lat_B_3_1,   Lat_B_3_2;
%         Lat_B_4_1;   Lat_B_4_2;
%        ];


%% Eigenvalues Computation

time            = (0 : 0.1 : 400)';
Long_StateVar    = zeros(4,length(time));
Lat_StateVar     = zeros(4,length(time));

%%%
% Longitudinal Motion Eigenvalue Computation
%%%

[Long_V,Long_D] = eig(Long_A);
    
    
%%%
% Lateral Motion Eigenvalue Computation
%%%

[Lat_V,Lat_D] = eig(Lat_A);


%% Plotting Transient Models

%%%
% Longitudinal Motion Transient Models
%%%

Long_lambda = diag(Long_D)
Long_C = Long_V\x_0;
% Long_C = inv(Long_V) * x_0;

for i = 1:length(time)
    Long_StateVar(:, i) =   Long_C(1) * Long_V(:,1) * exp(Long_lambda(1) * time(i)) + ...
                            Long_C(2) * Long_V(:,2) * exp(Long_lambda(2) * time(i)) + ...
                            Long_C(3) * Long_V(:,3) * exp(Long_lambda(3) * time(i)) + ...
                            Long_C(4) * Long_V(:,4) * exp(Long_lambda(4) * time(i));
end

figure(1);
plot(time, real(Long_StateVar(1,:)));
xlabel("Time (s)");
ylabel("\Deltau");

figure(2);
plot(time, real(Long_StateVar(2,:)));
xlabel("Time (s)");
ylabel("\Deltaw");

figure(3);
plot(time, real(Long_StateVar(3,:)));
xlabel("Time (s)");
ylabel("\Deltaq");

figure(4);
plot(time, real(Long_StateVar(4,:)));
xlabel("Time (s)");
ylabel("\Delta\theta");


%%%
% Lateral Motion Transient Models
%%%    

Lat_lambda   = diag(Lat_D);
Lat_C        = Lat_V\x_0;
% Lat_C        = inv(Lat_V) * x_0;

for i = 1:length(time)
    Lat_StateVar(:, i) =    Lat_C(1) * Lat_V(:,1) * exp(Lat_lambda(1) * time(i)) + ...
                            Lat_C(2) * Lat_V(:,2) * exp(Lat_lambda(2) * time(i)) + ...
                            Lat_C(3) * Lat_V(:,3) * exp(Lat_lambda(3) * time(i)) + ...
                            Lat_C(4) * Lat_V(:,4) * exp(Lat_lambda(4) * time(i));
end

figure(5);
plot(time, real(Lat_StateVar(1,:)));
xlabel("Time (s)");
ylabel("\Deltav");

figure(6);
plot(time, real(Lat_StateVar(2,:)));
xlabel("Time (s)");
ylabel("\Deltap");

figure(7);
plot(time, real(Lat_StateVar(3,:)));
xlabel("Time (s)");
ylabel("\Deltar");

figure(8);
plot(time, real(Lat_StateVar(4,:)));
xlabel("Time (s)");
ylabel("\Delta\phi");


%% Longitudinal Motion Solutions

%%%
% Finding Phugoid & Short Eigenvalues
%%%

if abs(real(Long_lambda(3))) > abs(real(Long_lambda(1)))
    
    Long_long_eta   = real(Long_lambda(1));
    Long_short_eta  = real(Long_lambda(3));
    
    Long_long_omega     = imag(Long_lambda(1));
    Long_short_omega    = imag(Long_lambda(3));

else
    
    Long_long_eta   = real(Long_lambda(3));
    Long_short_eta  = real(Long_lambda(1));
    
    Long_long_omega     = imag(Long_lambda(3));
    Long_short_omega    = imag(Long_lambda(1));

end

%%%
% Phugoid (Long) Period Solution
%%%

Long_long_t_half    = 0.69 / abs(Long_long_eta);
Long_long_period    = 2 * pi / Long_long_omega;
Long_long_N_half    = 0.110 * abs(Long_long_omega) / abs(Long_long_eta);

%%%
% Phugoid (Long) Period Solution Approximations
%%%

Long_omega_n_p          = sqrt(-Z_u * g / u_0);
Long_zeta_p             = - X_u / (2 * Long_omega_n_p);
Long_eigen_1_2_p_real   = - Long_zeta_p * Long_omega_n_p;
Long_eigen_1_2_p_imag   = Long_omega_n_p * sqrt(1 - Long_zeta_p ^ 2);
Long_period_p           = 2 * pi / Long_eigen_1_2_p_imag;
Long_t_half_p           = 0.69 / abs(Long_eigen_1_2_p_real);
Long_N_half_p           = 0.110 * Long_omega_n_p / abs(Long_eigen_1_2_p_real);

%%%
% Short Period Solution
%%%

Long_short_t_half       = 0.69 / abs(Long_short_eta);
Long_short_period       = 2 * pi / Long_short_omega;
Long_short_N_half       = 0.110 * abs(Long_short_omega) / abs(Long_short_eta);

%%%
% Short Period Solution Approximations
%%%

Long_omega_n_sp         = sqrt(Z_alpha * M_q / u_0 - M_alpha);
Long_zeta_sp            = - (M_q + M_alpha_dot + Z_alpha / u_0) / (2 * Long_omega_n_sp);
Long_eigen_1_2_sp_real  = - Long_zeta_sp * Long_omega_n_sp;
Long_eigen_1_2_sp_imag  = Long_omega_n_sp * sqrt(1 - Long_zeta_sp ^ 2);
Long_period_sp          = 2 * pi / Long_eigen_1_2_sp_imag;
Long_t_half_sp          = 0.69 / abs(Long_eigen_1_2_sp_real);
Long_N_half_sp          = 0.110 * Long_omega_n_sp / abs(Long_eigen_1_2_sp_real);

%% Table Generation

%%%
% Longitudinal Motion Approximation vs Exact Method Comparison Table
%%%

table_ApproxComp_col_Types = [""+...
    "Phugoid Half Amp Time (t_(1/2))";
    "Phugoid Period (P)";
    "Phugoid Frequency (omega)";
    "Short Half Amp Time (t_(1/2))";
    "Short Period (P)";
    "Short Frequency (omega)"];

table_ApproxComp_col_ExactMethod = [""+...
    num2str(Long_long_t_half) + " s";
    num2str(Long_long_period) + " s";
    num2str(Long_long_omega) + " rad/s"
    num2str(Long_short_t_half) + " s";
    num2str(Long_short_period) + " s";
    num2str(Long_short_omega) + " rad/s"];
           
table_ApproxComp_col_Approximation = [""+...
    num2str(Long_t_half_p) + " s";
    num2str(Long_period_p) + " s";
    num2str(Long_omega_n_p) + " rad/s"
    num2str(Long_t_half_sp) + " s";
    num2str(Long_period_sp) + " s";
    num2str(Long_omega_n_sp) + " rad/s"];
    
table_ApproxComp_col_Difference = [""+...
    num2str(100 * abs((Long_long_t_half - Long_t_half_p) / Long_long_t_half)) + "%";
    num2str(100 * abs((Long_long_period - Long_period_p) / Long_long_period)) + "%";
    num2str(100 * abs((Long_long_omega - Long_omega_n_p) / Long_long_omega)) + "%";
    num2str(100 * abs((Long_short_t_half - Long_t_half_sp) / Long_short_t_half)) + "%";
    num2str(100 * abs((Long_short_period - Long_period_sp) / Long_short_period)) + "%";
    num2str(100 * abs((Long_short_omega - Long_omega_n_sp) / Long_short_omega)) + "%"];

table_ApproxComp = table(table_ApproxComp_col_ExactMethod, table_ApproxComp_col_Approximation, table_ApproxComp_col_Difference);
table_ApproxComp.Properties.VariableNames = ["Exact Method", "Approximate Method", "Difference"];
table_ApproxComp.Properties.RowNames = table_ApproxComp_col_Types;

%%%
% Longitudinal Eigenvalue & Eigenvector Table
%%%

table_eigen_col = 0;



