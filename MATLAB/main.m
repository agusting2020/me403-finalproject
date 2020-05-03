
%
%   Control Systems Analysis of the
%   Transport Aircraft: Convair 880
%   
%
%   Created by: Agustin Guzman
%

%% Import Aircraft

% aircraft_sample

aircraft_convair_880

%% Missing Lateral Motion Coefficients

C_y_p       = 0;
C_y_r       = 0;
C_y_delta_a = 0;    


%% Initializing Lateral Motion Derivatives

%%%
% Directional Derivatives
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

L_v         = L_beta / u_0;
N_v         = N_beta / u_0;
Y_v         = Y_beta / u_0;

%%%
% Star Directional Derivatives
%%%

L_star_v    = L_v / (1 - (I_xz ^ 2 / (I_x * I_z)));
N_star_v    = N_v / (1 - (I_xz ^ 2 / (I_x * I_z)));

L_star_r    = L_r / (1 - (I_xz ^ 2 / (I_x * I_z)));
N_star_r    = N_r / (1 - (I_xz ^ 2 / (I_x * I_z)));

L_star_p    = L_p / (1 - (I_xz ^ 2 / (I_x * I_z)));
N_star_p    = N_p / (1 - (I_xz ^ 2 / (I_x * I_z)));

L_star_delta_a   = L_delta_a / (1 - (I_xz ^ 2 / (I_x * I_z)));
N_star_delta_a   = N_delta_a / (1 - (I_xz ^ 2 / (I_x * I_z)));

N_star_delta_r  = N_delta_r / (1 - (I_xz ^ 2 / (I_x * I_z)));
L_star_delta_r  = L_delta_r / (1 - (I_xz ^ 2 / (I_x * I_z)));


%% Calculating Lateral Motion Matrices

%%%
% Lateral Motion Matrix A
%%%

A_1_1 = Y_v;
A_1_2 = Y_p;
A_1_3 = -(u_0 - Y_r);
A_1_4 = g * cos(theta_0);

A_2_1 = L_star_v + (I_xz / I_x) * N_star_v;
A_2_2 = L_star_p + (I_xz / I_x) * N_star_p;
A_2_3 = L_star_r + (I_xz / I_x) * N_star_r;
A_2_4 = 0;

A_3_1 = N_star_v + (I_xz / I_z) * L_star_v;
A_3_2 = N_star_p + (I_xz / I_z) * L_star_p;
A_3_3 = N_star_r + (I_xz / I_z) * L_star_r;
A_3_4 = 0;

A_4_1 = 0;
A_4_2 = 1;
A_4_3 = 0;
A_4_4 = 0;

A = [
    A_1_1,   A_1_2,   A_1_3,   A_1_4;
    A_2_1,   A_2_2,   A_2_3,   A_2_4;
    A_3_1,   A_3_2,   A_3_3,   A_3_4;
    A_4_1,   A_4_2,   A_4_3,   A_4_4]


%%%
% Lateral Motion Matrix B
%%%

B_1_1 = 0;
B_1_2 = Y_delta_r;

B_2_1 = L_star_delta_a + (I_xz / I_x) * N_star_delta_a;
B_2_2 = L_star_delta_r + (I_xz / I_x) * N_star_delta_r;

B_3_1 = N_star_delta_a + (I_xz / I_x) * L_star_delta_a;
B_3_2 = N_star_delta_r + (I_xz / I_x) * L_star_delta_r;

B_4_1 = 0;
B_4_2 = 0;

B = [
    B_1_1,   B_1_2;
    B_2_1,   B_2_2;
    B_3_1,   B_3_2;
    B_4_1,   B_4_2]


%% Computing Lateral Motion Eigenvalues

[V,D]   = eig(A)
lambda  = diag(D)

eta_roll    = real(lambda(1));
eta_dutch   = real(lambda(2));
eta_spiral  = real(lambda(4));

omega_dutch_pos  = imag(lambda(2));
omega_dutch_neg  = imag(lambda(3));

%% Calculating Lateral Motion Approximations

%%%
% Spiral Approximation
%%%

approx_lambda_spiral    = (L_beta * N_r - L_r * N_beta) / L_beta;
approx_eta_spiral       = real(approx_lambda_spiral);

%%%
% Roll Approximation
%%%

approx_lambda_roll = L_p;
approx_eta_roll       = real(approx_lambda_roll);

%%%
% Dutch Approximation
%%%

approx_omega_n_dutch = sqrt((Y_beta * N_r - N_beta * Y_r + u_0 * N_p) / u_0);

approx_zeta_dutch = - (1 / (2 * approx_omega_n_dutch)) * ((Y_beta + u_0 * N_r) / u_0);





%% Plotting Lateral Motion Transient Models

time        = (0 : 0.1 : 100)';
StateVar    = zeros(4,length(time));
C           = V\x_0;

for i = 1:length(time)
    StateVar(:, i)  =   C(1) * V(:,1) * exp(lambda(1) * time(i)) + ...
                        C(2) * V(:,2) * exp(lambda(2) * time(i)) + ...
                        C(3) * V(:,3) * exp(lambda(3) * time(i)) + ...
                        C(4) * V(:,4) * exp(lambda(4) * time(i));
end

figure(1);
plot(time, real(StateVar(1,:)));
xlabel("Time (s)");
ylabel("\Deltav");

figure(2);
plot(time, real(StateVar(2,:)));
xlabel("Time (s)");
ylabel("\Deltap");

figure(3);
plot(time, real(StateVar(3,:)));
xlabel("Time (s)");
ylabel("\Deltar");

figure(4);
plot(time, real(StateVar(4,:)));
xlabel("Time (s)");
ylabel("\Delta\phi");


%% Lateral Motion Table Generation

%%%
% Lateral Motion Approximation vs Exact Method Comparison Table
%%%

table_ApproxComp_col_Types = [""+...
    "Phugoid Half Amp Time (t_(1/2))";
    "Phugoid Period (P)";
    "Phugoid Frequency (omega)";
    "Short Half Amp Time (t_(1/2))";
    "Short Period (P)";
    "Short Frequency (omega)"];

table_ApproxComp_col_ExactMethod = [""+...
    num2str(long_t_half) + " s";
    num2str(long_period) + " s";
    num2str(long_omega) + " rad/s"
    num2str(short_t_half) + " s";
    num2str(short_period) + " s";
    num2str(short_omega) + " rad/s"];
           
table_ApproxComp_col_Approximation = [""+...
    num2str(t_half_p) + " s";
    num2str(period_p) + " s";
    num2str(omega_n_p) + " rad/s"
    num2str(t_half_sp) + " s";
    num2str(period_sp) + " s";
    num2str(omega_n_sp) + " rad/s"];
    
table_ApproxComp_col_Difference = [""+...
    num2str(100 * abs((long_t_half - t_half_p) / long_t_half)) + "%";
    num2str(100 * abs((long_period - period_p) / long_period)) + "%";
    num2str(100 * abs((long_omega - omega_n_p) / long_omega)) + "%";
    num2str(100 * abs((short_t_half - t_half_sp) / short_t_half)) + "%";
    num2str(100 * abs((short_period - period_sp) / short_period)) + "%";
    num2str(100 * abs((short_omega - omega_n_sp) / short_omega)) + "%"];

table_ApproxComp = table(table_ApproxComp_col_ExactMethod, table_ApproxComp_col_Approximation, table_ApproxComp_col_Difference);
table_ApproxComp.Properties.VariableNames = ["Exact Method", "Approximate Method", "Difference"];
table_ApproxComp.Properties.RowNames = table_ApproxComp_col_Types;


%% Part II Transfer Function

%%% 
% Computing Transfer Function G(s)
%%%

compute_G

%%%
% Obtaining Transfer Function Information
%%%

sysinfo(G)

%%%
% Plotting Step Response
%%%

figure(9);
step(G)

%%%
% Setting Up Proportional Controller
%%%

k_p = 0.1;
G_c = k_p;

% * (1 + (1 / T_i * s) + T_d * s);

%%%
% Plotting Root Locus
%%%

figure(10);
rlocus(G_c * G)

%%%
% 


