
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
    A_4_1,   A_4_2,   A_4_3,   A_4_4];


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
    B_4_1,   B_4_2];


%% Computing Lateral Motion Eigenvalues

[V,D]   = eig(A);
lambda  = diag(D);

% eta_roll    = real(lambda(1));
% eta_dutch   = real(lambda(2));
% eta_spiral  = real(lambda(4));
% 
% omega_dutch_pos  = imag(lambda(2));
% omega_dutch_neg  = imag(lambda(3));

lambda_roll = lambda(1);
lambda_dutch = lambda(2);
lambda_dutch_negative = lambda(3);
lambda_spiral = lambda(4);

eta_roll    = real(lambda_roll);
eta_dutch   = real(lambda_dutch);
eta_spiral  = real(lambda_spiral);

omega_dutch_pos  = imag(lambda_dutch);
omega_dutch_neg  = imag(lambda_dutch_negative);


%% Computing Lateral Motion Solutions

period_dutch    = 2 * pi / abs(omega_dutch_pos);

tau_roll   = - 1 / eta_roll;
tau_dutch  = - 1 / eta_dutch;
tau_spiral = - 1 / eta_spiral;

time_half_roll     = 0.69 / abs(eta_roll);
time_half_dutch    = 0.69 / abs(eta_dutch);
time_half_spiral   = 0.69 / abs(eta_spiral);



%% Computing Lateral Motion Eigenvalue & Solution Approximations

%%%
% Spiral Eigenvalue Approximation
%%%

approx_lambda_spiral    = (L_beta * N_r - L_r * N_beta) / L_beta;
approx_eta_spiral       = real(approx_lambda_spiral);


%%%
% Roll Eigenvalue Approximation
%%%

approx_lambda_roll  = L_p;
approx_eta_roll     = real(approx_lambda_roll);

%%%
% Dutch Eigenvalue Approximation
%%%

eqa = 1;
eqb = -(Y_beta + u_0 * N_r) / u_0;
eqc = (Y_beta * N_r - N_beta * Y_r + u_0 * N_beta) / u_0;

approx_lambda_dutch = - eqb / (2 * eqa) - sqrt((eqb/ (2 * eqa)) ^ 2 - eqc / eqa);

% approx_omega_n_dutch = sqrt((Y_beta * N_r - N_beta * Y_r + u_0 * N_beta) / u_0);
% approx_zeta_dutch    = -(0.5 / approx_omega_n_dutch) * ((Y_beta + u_0 * N_r) / u_0);

approx_omega_dutch = abs(imag(approx_lambda_dutch));
approx_eta_dutch = abs(real(approx_lambda_dutch));

approx_omega_n_dutch = sqrt((Y_beta * N_r - N_beta * Y_r + u_0 * N_beta) / u_0);
approx_zeta_dutch    = -(0.5 / approx_omega_dutch) * ((Y_beta + u_0 * N_r) / u_0);

approx_eta_dutch_calc = -approx_omega_n_dutch * approx_zeta_dutch;
approx_omega_dutch_calc = approx_omega_n_dutch * sqrt(1 - approx_zeta_dutch ^ 2);


%% Lateral Motion Solution Approximations

approx_period_dutch    = 2 * pi / abs(approx_omega_dutch);

approx_tau_roll   = - 1 / approx_eta_roll;
approx_tau_dutch  = - 1 / approx_eta_dutch;
approx_tau_spiral = - 1 / approx_eta_spiral;

approx_time_half_roll     = 0.69 / abs(approx_eta_roll);
approx_time_half_dutch    = 0.69 / abs(approx_eta_dutch);
approx_time_half_spiral   = 0.69 / abs(approx_eta_spiral);


%% Plotting Lateral Motion Transient Models

time        = (0 : 0.1 : 100)';
StateVar    = zeros(4,length(time));
C           = V\x_0;
StateVar2   = StateVar;
C2          = V\[0; 1; 0; (15*pi/180)];

for i = 1:length(time)
    StateVar(:, i)  =   C(1) * V(:,1) * exp(lambda(1) * time(i)) + ...
                        C(2) * V(:,2) * exp(lambda(2) * time(i)) + ...
                        C(3) * V(:,3) * exp(lambda(3) * time(i)) + ...
                        C(4) * V(:,4) * exp(lambda(4) * time(i));
end

for i = 1:length(time)
    StateVar2(:, i) =   C2(1) * V(:,1) * exp(lambda(1) * time(i)) + ...
                        C2(2) * V(:,2) * exp(lambda(2) * time(i)) + ...
                        C2(3) * V(:,3) * exp(lambda(3) * time(i)) + ...
                        C2(4) * V(:,4) * exp(lambda(4) * time(i));
end

figure(1);

    tiledlayout(2,2)

        nexttile
        plot(time, real(StateVar(1,:)));
        title("Transient Response for \Deltav")
        xlabel("Time (s)");
        ylabel("\Deltav");

        nexttile
        plot(time, real(StateVar(2,:)));
        title("Transient Response for \Deltap")
        xlabel("Time (s)");
        ylabel("\Deltap");

        nexttile
        plot(time, real(StateVar(3,:)));
        title("Transient Response for \Deltar")
        xlabel("Time (s)");
        ylabel("\Deltar");

        nexttile
        plot(time, real(StateVar(4,:)));
        title("Transient Response for \Delta\phi")
        xlabel("Time (s)");
        ylabel("\Delta\phi");
        
        
figure(2);

    tiledlayout(2,2)

        nexttile
        plot(time, real(StateVar2(1,:)));
        title("Transient Response for \Deltav")
        xlabel("Time (s)");
        ylabel("\Deltav");

        nexttile
        plot(time, real(StateVar2(2,:)));
        title("Transient Response for \Deltap")
        xlabel("Time (s)");
        ylabel("\Deltap");

        nexttile
        plot(time, real(StateVar2(3,:)));
        title("Transient Response for \Deltar")
        xlabel("Time (s)");
        ylabel("\Deltar");

        nexttile
        plot(time, real(StateVar2(4,:)));
        title("Transient Response for \Delta\phi")
        xlabel("Time (s)");
        ylabel("\Delta\phi");


%% Lateral Motion Table Generation

%%%
% Eigenvectors
%%%

table_Eigenvectors_col_Types = [""+...
    "\Deltav Eigenvector";
    "\Deltap Eigenvector";
    "\Deltar Eigenvector"];

table_Eigenvectors_col_Eigenvector_Roll = [""+...
    num2str(V(1,1));
    num2str(V(2,1));
    num2str(V(3,1))];

table_Eigenvectors_col_Eigenvector_Dutch_1 = [""+...
    num2str(V(1,2));
    num2str(V(2,2));
    num2str(V(3,2))];

table_Eigenvectors_col_Eigenvector_Dutch_2 = [""+...
    num2str(V(1,3));
    num2str(V(2,3));
    num2str(V(3,3))];

table_Eigenvectors_col_Eigenvector_Spiral = [""+...
    num2str(V(1,4));
    num2str(V(2,4));
    num2str(V(3,4))];

table_Eigenvectors = table(table_Eigenvectors_col_Eigenvector_Roll, table_Eigenvectors_col_Eigenvector_Dutch_1, table_Eigenvectors_col_Eigenvector_Dutch_2, table_Eigenvectors_col_Eigenvector_Spiral);
table_Eigenvectors.Properties.VariableNames = ["Roll Eigenvector", "Dutch Roll Eigenvector 1", "Dutch Roll Eigenvector 2", "Spiral Eigenvector"];
table_Eigenvectors.Properties.RowNames = table_Eigenvectors_col_Types;

%%% print table
writetable(table_Eigenvectors,'tabulation.xlsx','Sheet',1);
table_Eigenvectors


%%%
% Lateral Motion Approximation vs Exact Method Comparison Table
%%%

table_ApproxComp_col_Types = [""+...
    "Roll Eigenvalue";
    "Dutch Roll Eigenvalue";
    "Spiral Eigenvalue";
    "Roll Time Constant (s)";
    "Roll Half Amp Time (t_(1/2))";
    "Roll Period (P)";
    "Roll Frequency (omega)";
    "Dutch Roll Time Constant (s)";
    "Dutch Roll Half Amp Time (t_(1/2))";
    "Dutch Roll Period (P)";
    "Dutch Roll Frequency (omega)";
    "Spiral Time Constant (s)";
    "Spiral Half Amp Time (t_(1/2))";
    "Spiral Period (P)";
    "Spiral Frequency (omega)"];

table_ApproxComp_col_ExactMethod = [""+...
    num2str(lambda_roll);
    num2str(lambda_dutch_negative);
    num2str(lambda_spiral);
    num2str(tau_roll) + " s";
    num2str(time_half_roll) + " s";
    "-"; % no period for roll
    "-"; % no frequency for roll
    num2str(tau_dutch) + " s";
    num2str(time_half_dutch) + " s";
    num2str(period_dutch) + " s";
    num2str(omega_dutch_pos) + " rad/s"
    num2str(tau_spiral) + " s";
    num2str(time_half_spiral) + " s";
    "-"; % no period for spiral
    "-"]; % no frequency for spiral
           
table_ApproxComp_col_Approximation = [""+...
    num2str(approx_lambda_roll);
    num2str(approx_lambda_dutch);
    num2str(approx_lambda_spiral);
    
    num2str(approx_tau_roll) + " s";
    num2str(approx_time_half_roll) + " s";
    "-"; % no period for roll
    "-"; % no frequency for roll
    
    num2str(approx_tau_dutch) + " s";
    num2str(approx_time_half_dutch) + " s";
    num2str(approx_period_dutch) + " s";
    num2str(approx_omega_dutch) + " rad/s"
    
    num2str(approx_tau_spiral) + " s";
    num2str(approx_time_half_spiral) + " s";
    "-"; % no period for roll
    "-"]; % no frequency for roll
    
table_ApproxComp_col_Difference = [""+...
    num2str(100 * abs((lambda_roll - approx_lambda_roll) / lambda_roll)) + "%";
    num2str(100 * abs((lambda_dutch - approx_lambda_dutch) / lambda_dutch)) + "%";
    num2str(100 * abs((lambda_spiral - approx_lambda_spiral) / lambda_spiral)) + "%";
    
    num2str(100 * abs((tau_roll - approx_tau_roll) / tau_roll)) + "%";
    num2str(100 * abs((time_half_roll - approx_time_half_roll) / time_half_roll)) + "%";
    "-"; % no period for roll
    "-"; % no frequency for roll
    
    num2str(100 * abs((tau_dutch - approx_tau_dutch) / tau_dutch)) + "%";
    num2str(100 * abs((time_half_dutch - approx_time_half_dutch) / time_half_dutch)) + "%";
    num2str(100 * abs((period_dutch - approx_period_dutch) / period_dutch)) + "%";
    num2str(100 * abs((omega_dutch_pos - approx_omega_dutch) / omega_dutch_pos)) + "%";
    
    num2str(100 * abs((tau_spiral - approx_tau_spiral) / tau_spiral)) + "%";
    num2str(100 * abs((time_half_spiral - approx_time_half_spiral) / time_half_spiral)) + "%";
    "-"; % no period for spiral
    "-"]; % no frequency for spiral

table_ApproxComp = table(table_ApproxComp_col_ExactMethod, table_ApproxComp_col_Approximation, table_ApproxComp_col_Difference);
table_ApproxComp.Properties.VariableNames = ["Exact Method", "Approximate Method", "Difference"];
table_ApproxComp.Properties.RowNames = table_ApproxComp_col_Types;

%%% print table
writetable(table_ApproxComp,'tabulation.xlsx','Sheet',2);
table_ApproxComp


%% Part II Transfer Function

%%% 
% Computing Transfer Function G(s)
%%%

compute_G


%%% print transient function G
G


%%%
% Obtaining Transfer Function Information
%%%

stepinfo(G)


%%%
% Plotting Step Responses & Root Locus
%%%

responses_aircraft = figure(3);

    tiledlayout(2,2)
        
        nexttile([1 2])
        plot(step(G))
        title("Step Response Without Controller")
        xlabel("Time (s)");
        ylabel("Amplitude");        

        nexttile([1 2])
        plot(impulse(G))
        title("Impulse Response Without Controller")
        xlabel("Time (s)");
        ylabel("Amplitude");


%% Proportional Controller

%%%
% Defining Parameters
%%%

% C_P_PT = pidtune(G,'p');

% [Kcr_P,r_P] = rlocfind(G)
% 
% Kp_P = Kcr_P * 0.5

Kp_P = [-0.385429347282575];

Kp_p = 0.05

C_P = pid(Kp_P);

T_P = feedback(G*C_P,1);


%%%
% Plotting Proportional Controller Response & Root Locus
%%%

responses_controller = figure(11);

    tiledlayout(2,2)
        
        nexttile([1 2])
        plot(step(T_P))
        title("Step Response With Proportional Controller")
        xlabel("Time (s)");
        ylabel("Amplitude");        

        nexttile([1 2])
        plot(impulse(T_P))
        title("Impulse Response With Proportional Controller")
        xlabel("Time (s)");
        ylabel("Amplitude");

root_locus = figure(12);

    rlocus(C_P * G)
    title("Root Locus (Proportional Controller)");
%     sgrid(sqrt(2)/2,)
    axis equal


%% PID Controller

%%%
% Defining Parameters for PID Controller
%%%

% Kp_PI = -0.00018;
% 
% Ki_PI = -0.0000556;
% 
% C_PI  = pid(Kp_PI, Ki_PI);

C_PID = pidtune(G,'pid')

T_PID = feedback(G*C_PID,1); 

% * (1 + (1 / T_i * s) + T_d * s);


%%%
% Plotting PID Controller  Response & Root Locus
%%%
    
responses_controller_PID = figure(21);

    tiledlayout(2,2)
        
        nexttile([1 2])
        plot(step(T_PID))
        title("Step Response With PID Controller")
        xlabel("Time (s)");
        ylabel("Amplitude");        

        nexttile([1 2])
        plot(impulse(T_PID))
        title("Impulse Response With PID Controller")
        xlabel("Time (s)");
        ylabel("Amplitude");

root_locus_PID = figure(22);

    rlocus(T_PID)
    title("Root Locus (PID Controller)");
    
    
% [Kcr, r] = rlocfind(T_PID)
% 
% Tcr = 2*pi/imag(r(2));
% 
% Kp=0.6*Kcr;
% 
% Ti=0.5*Tcr;
% 
% Td=0.125*Tcr;
% 
% C_PID = tf([Kp*Td*Ti Kp*Ti Kp], [Ti 0])
% 
% T_PID = feedback(G*C_PID,1);
% 
% root_locus_PID_2 = figure(23);
% 
%     rlocus(T_PID)
%     title("Root Locus (PID Controller)");



