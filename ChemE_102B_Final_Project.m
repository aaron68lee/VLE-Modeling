%{%% %%%%%%%%%%%%%%%%%%%% Header %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Name: Aaron Lee
UID: 505 540 473
Date: 6/10/22
Objective: Simulate Vapor-Liquid Equilibrium (VLE) Systems

Water BubbleT = 81.6886;                liquid  vapor     liquid vapor
Water: 0.1052, 0.2079, 0.1679, 0.0776, [0.1387, 0.2428], [0.21, 0.1]
%}

clc; clear all;

%% %%%%%%%%%%%%%%%%%%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%

MAX = 4; % number of species

SOL_q2a_T = -1; % bubble point temp, use iteration 
SOL_q2a_y = zeros(1, MAX); % bubble composition array (1x4)

SOL_q2b_T = -1; % dew point temp, use iteration
SOL_q2b_x = zeros(1, MAX); % dew composition array (1x4)

% use virial equation of state and UNIFAC
SOL_q3_bublT = -1; % bubble point temp, use iteration
SOL_q3_dewT = -1; % dew point temp
SOL_q3_y = zeros(1, MAX); % bubble composition
SOL_q3_x = zeros(1, MAX); % dew composition

% VLE compositions
SOL_q4_y = zeros(1, MAX); % vapor phase composition
SOL_q4_x = zeros(1, MAX); % liquid phase composition
SOL_q4_V = -1; % vapor mole fraction
SOL_q4_L = -1; % vapor mole fraction

% nonideal mixtures of solutions and gases with fugacity coefficients
SOL_q5_y = zeros(1, MAX); % vapor phase composition
SOL_q5_T = -1; % temp
SOL_q5_P = -1; % pressure
SOL_q5_V = -1;
SOL_q5_L = -1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ERROR = 1e-9; % max tolerance error limit
R = 8.3145; % J/mol*K
comp = [0.21, 0.19, 0.27, 0.33];
SIZE = size(comp, 2);

antoine = [16.3872, 3885.70, 230.170; %or use space and "..." after line
           16.8958, 3795.17, 230.918;
           14.3145, 2756.22, 228.060;
           15.3144, 3213.43, 182.739];
       
A = antoine(:, 1);
B = antoine(:, 2);
C = antoine(:, 3);

%% %%%%%%%%%%%%%%%%%%%%% Debugging Section %%%%%%%%%%%%%%%%%%%%%%%

Psat = zeros(1, MAX);
for i = 1:MAX
    % find saturated pressure
    Psat(i) = pressure(A(i), B(i), C(i), 100);
end

phiTest = Phi([0.1, 0.2, 0.3, 0.4], antoine, 101.325, 298.15);
%disp(phiTest);

gammaTest = activity([0.1, 0.2, 0.3, 0.4], 298);
%disp(gammaTest);

%% %%%%%%%%%%%%%%%%%%%%% Question 2a Bubble Temp %%%%%%%%%%%%%%%%%%%%%%%

Pre = 101.325; % pressure in kPa
P = zeros(1, MAX);
temp = -1;

Tsat = B ./ (A - log(Pre)) - C; % 1 x MAX array
SOL_q2a_T = sum(comp * Tsat);

while abs(temp - SOL_q2a_T) > ERROR
    
    temp = SOL_q2a_T;
    for i = 1:MAX
        % find saturated pressure
        P(i) = pressure(A(i), B(i), C(i), SOL_q2a_T);
    end

    Pjsat = Pre / (sum(comp .* P / P(1)));
    
    SOL_q2a_T = temperature(A(1), B(1), C(1), Pjsat);
    SOL_q2a_y = comp .* P / Pre;

end

fprintf("Q2a Bubble Temp: " + SOL_q2a_T + "\n");
fprintf("Q2a Bubble Point Comp: ");
disp(SOL_q2a_y);

%% %%%%%%%%%%%%%%%%%%%%% Question 2b %%%%%%%%%%%%%%%%%%%%%%%

Pre = 101.325; % pressure in kPa
P = zeros(1, MAX);
temp = -1;

Tsat = B ./ (A - log(Pre)) - C;
SOL_q2b_T = sum(comp * Tsat); % initial guess value for temp
for i = 1:MAX
    P(i) = pressure(A(i), B(i), C(i), SOL_q2b_T);
end
Pjsat = Pre * (sum(comp .* P(1) ./ P)); % comp in gas phase: y
SOL_q2b_T = temperature(A(1), B(1), C(1), Pjsat);

SOL_q2b_x = comp .* Pre ./ (P); % comp in liquid: x

Pjsat = Pre * (sum(comp .* P(1) ./ P)); % comp in gas phase: y

while abs(temp - SOL_q2b_T) > ERROR
    
    for i = 1:MAX
        % find saturated pressure
        P(i) = pressure(A(i), B(i), C(i), SOL_q2b_T);
    end
    
    % calculate dew point comp   
    SOL_q2b_x = comp .* Pre ./ (P); % comp in liquid: x
    SOL_q2b_x = SOL_q2b_x ./ (sum(SOL_q2b_x));
    
    temp = SOL_q2b_T;
    Pjsat = Pre * (sum(comp .* P(1) ./ P)); % comp in gas phase: y
    SOL_q2b_T = temperature(A(1), B(1), C(1), Pjsat);
    
end

fprintf("Q2b Dew Temp: " + SOL_q2b_T + "\n");
fprintf("Q2b Dew Point Comp: ");
disp(SOL_q2b_x);

%% %%%%%%%%%%%%%%%%%%%%% Question 3 %%%%%%%%%%%%%%%%%%%%%%%

%% find bubble point temperature & comp NONIDEAL

temp = -1;
phi = ones(1, MAX); % initial guess

Tsat = B ./ (A - log(Pre)) - C;
SOL_q3_bublT = sum(comp * Tsat); % initial guess value for temp
% find saturated pressure
for i = 1:MAX
    P(i) = pressure(A(i), B(i), C(i), SOL_q3_bublT);
end
activities = activity(comp, SOL_q3_bublT);
Pjsat = Pre / (sum(comp .* activities ./ phi .* P ./ P(1)));
SOL_q3_bublT = temperature(A(1), B(1), C(1), Pjsat);

% iterative solver
while abs(temp - SOL_q3_bublT) > ERROR
    
    % find saturated pressure
    for i = 1:MAX
        P(i) = pressure(A(i), B(i), C(i), SOL_q3_bublT);
    end
    
    % update variable guesses
    SOL_q3_y = comp .* activities .* P ./ (phi .* Pre); % comp in liquid: x
    
    % propogation
    %{SOL_q3_y = [0.1895    0.1787    0.6698    0.0630];
    %T = 74.7545 + 273.15;
    phi = Phi(SOL_q3_y, antoine, Pre, SOL_q3_bublT + 273.15); 
    activities = activity(comp, SOL_q3_bublT);
    Pjsat = Pre / (sum(comp .* activities ./ phi .* P ./ P(1)));
    
    temp = SOL_q3_bublT;
    SOL_q3_bublT = temperature(A(1), B(1), C(1), Pjsat);
    
end

fprintf("Q3 Bubble Temp: " + SOL_q3_bublT + "\n");
fprintf("Q3 Bubble Comp: ");
disp(SOL_q3_y);

%% find dew point temperature & comp NONIDEAL

temp = -1;
temp_gamma = -1 * ones(1, MAX);
phi = ones(1, MAX); % initial guess
activities = ones(1, MAX);

Tsat = B ./ (A - log(Pre)) - C;
SOL_q3_dewT = sum(comp * Tsat); % initial guess value for temp
for i = 1:MAX
    P(i) = pressure(A(i), B(i), C(i), SOL_q3_dewT);
end
Pjsat = Pre * (sum(P(1) .* comp ./ P)); % comp in gas phase: y
SOL_q3_dewT = temperature(A(1), B(1), C(1), Pjsat);
for i = 1:MAX
    P(i) = pressure(A(i), B(i), C(i), SOL_q3_dewT);
end
phi = Phi(comp, antoine, Pre, SOL_q3_bublT + 273.15);
SOL_q3_x = comp .* phi * Pre ./ (activities .* P); % comp in liquid: x
activities = activity(SOL_q3_x, SOL_q3_dewT);

Pjsat = Pre * (sum(comp .* phi ./ activities .* P(1) ./ P)); % comp in gas phase: y
SOL_q3_dewT = temperature(A(1), B(1), C(1), Pjsat);

% iterative solver
while abs(temp - SOL_q3_dewT) > ERROR

    % find saturated pressure
    for i = 1:MAX
        P(i) = pressure(A(i), B(i), C(i), SOL_q3_dewT);
    end
    
    % update variable guesses
    phi = Phi(SOL_q3_y, antoine, Pre, SOL_q3_dewT + 273.15);
    
    % %%%%%%%%%%% propogation %%%%%%%%%%%
    
    % normalization
    while ~gammaCondition(activities, temp_gamma)
        temp_gamma = activities;
        
        SOL_q3_x = comp .* phi * Pre ./ (activities .* P); % comp in liquid: x
        % normalize x vector
        SOL_q3_x = SOL_q3_x ./ (sum(SOL_q3_x));
        %SOL_q3_x2 = SOL_q3_x / (norm(SOL_q3_x));
        %x = x ./ (sum(x));
        activities = activity(SOL_q3_x, SOL_q3_dewT);
    end
    
    temp = SOL_q3_dewT;
    Pjsat = Pre * (sum(comp .* phi ./ activities .* P(1) ./ P)); % comp in gas phase: y
    SOL_q3_dewT = temperature(A(1), B(1), C(1), Pjsat);
    
end

fprintf("Q3 Dew Temp: " + SOL_q3_dewT + "\n");
fprintf("Q3 Dew Comp: ");
disp(SOL_q3_x);

%% %%%%%%%%%%%%%%%%%%%%% Question 4 %%%%%%%%%%%%%%%%%%%%%%%

T = 360; % in Kelvin
Pre = 101.325; % in kPa
    
%% Bubble and Dew Composition Sub-calculations

%% Dew Point Pressure

temp = -1;
temp_gamma = -1 * ones(1, MAX);
x = zeros(1, MAX);
phi = ones(1, MAX); % initial guess
activities = ones(1, MAX);
% find saturated pressure
for i = 1:MAX
    P(i) = pressure(A(i), B(i), C(i), T - 273.15);
end
dewP = 1 / sum(comp .* phi ./ (activities .* P));
x = comp .* phi .* Pre ./ (activities .* P);
activities = activity(x, T - 273.15);
dewP = 1 / sum(comp .* phi ./ (activities .* P));

% iterative dewP solver
while abs(temp - dewP) > ERROR
    
    phi = Phi(comp, antoine, dewP, T); % pass T in Kelvin 
    
    % normalization
    while ~gammaCondition(activities, temp_gamma) % check each gamma < ERROR
        temp_gamma = activities;
        
        x = comp .* phi .* Pre ./ (activities .* P); % comp in liquid: x
        % normalize x vector
        x = x ./ (sum(x));
        %x = x / (norm(x));
        activities = activity(x, T - 273.15); % pass T in Celsius 
    end
    temp = dewP;
    dewP = 1 / sum(comp .* phi ./ (activities .* P));
end

%% Bubble Point Pressure

temp = -1;
y = zeros(1, MAX);
phi = ones(1, MAX); % initial guess
activities = activity(comp, T - 273.15);
% find saturated pressure
for i = 1:MAX
    P(i) = pressure(A(i), B(i), C(i), T - 273.15);
end
bubP = sum(comp .* activities .* P ./ phi);

% iterative bubP solver
while abs(temp - bubP) > ERROR
    y = comp .* activities .* P ./ (phi .* bubP);
    phi = Phi(y, antoine, bubP, T);
    temp = bubP;
    bubP = sum(comp .* activities .* P ./ phi);
end

%% VLE Solver
tempX = -1 * ones(1, MAX);
tempY = -1 * ones(1, MAX);
tempV = -1;
V1 = -1;
MAX_ITER = 1000;

if (Pre > dewP) && (Pre < bubP)
    V1 = 0.5;
    Varray(1) = V1; % initial point for V
    
    while ~VLE_Condition(abs(tempV-Varray(end)), abs(tempX-x), abs(tempY-y))
        K1 = activities .* P ./ (phi .* Pre); %y ./ x; % alternative K = ;
        z1 = comp;
        F = 0.5;
        % Newton's method for dF/dV
        %syms f(z, K, V)
        %syms g(z, K, V)
        f = @(z, K, V) sum((z .* (K - 1)) ./ (1 + V .* (K - 1)));
        g = @(z, K, V) -1*sum( (z .* (K - 1).^2) ./ (1 + V .* (K - 1) ).^2 ); % g = df/dV
        
        counter = 1;
        tempV = Varray(end);
        % Newton's Method
        while abs(F) > ERROR
            counter = counter + 1;
            F = f(z1, K1, Varray(end));
            dF = g(z1, K1, Varray(end));
            Varray(end+1) = Varray(end) - F / dF;
            if (counter > MAX_ITER)
                fprintf("\nMax iterations reached without convergence.\n");
                break;
            end
        end

        tempX = x;
        tempY = y;
        
        x = z1 ./ (1 + Varray(end) .* (K1 - 1));
        y = K1 .* x;
        activities = activity(x, T - 273.15); 
        phi = Phi(y, antoine, Pre, T); % set breakpoint here
    end
    
    %% Print final compositions (dew, bubble)
    
    SOL_q4_y = y;
    SOL_q4_x = x;
    SOL_q4_V = Varray(end);
    SOL_q4_L = 1 - SOL_q4_V;
    
    %% find V-L comp

    fprintf("Q4 Vapor Fraction: " + SOL_q4_V + "\n");
    fprintf("Q4 Liquid Fraction: " + SOL_q4_L + "\n");
    
    fprintf("Q4 Dew Liquid Comp: ");
    disp(SOL_q4_x);
    fprintf("Q4 Bubble Vapor Comp: ");
    disp(SOL_q4_y);

else
    fprintf("\nVLE Requirement Not Met\n");
end

%% %%%%%%%%%%%%%%%%%%%%% Question 5 %%%%%%%%%%%%%%%%%%%%%%%

%% find equilibrium point temperature & pressure & comp NONIDEAL

temp = -1;
phi = ones(1, MAX); % initial guess

Tsat = B ./ (A - log(Pre)) - C;
SOL_q5_T = sum(comp * Tsat); % initial guess value for temp
% find saturated pressure
for i = 1:MAX
    P(i) = pressure(A(i), B(i), C(i), SOL_q5_T);
end
activities = activity(comp, SOL_q5_T);
Pjsat = Pre / (sum(comp .* activities ./ phi .* P ./ P(1)));
SOL_q5_T = temperature(A(1), B(1), C(1), Pjsat);
SOL_q5_y = [0.10, 0.20, 0.30, 0.40]; % random guess with 0.1 water given

% iterative solver
while abs(temp - SOL_q5_T) > ERROR
    
    % find saturated pressure
    for i = 1:MAX
        P(i) = pressure(A(i), B(i), C(i), SOL_q5_T);
    end
    
    SOL_q5_P = sum(P .* comp .* activities);
    % update variable guesses
    SOL_q5_y = comp .* activities .* P ./ (phi .* SOL_q5_P); % comp in liquid: x
    %SOL_q5_y(1) = 0.1;
    
    % propogation
    phi = Phi(SOL_q5_y, antoine, SOL_q5_P, SOL_q5_T + 273.15); %T;
    activities = activity(comp, SOL_q2a_T);
    Pjsat = SOL_q5_P / (sum(comp .* activities ./ phi .* P ./ P(1)));
    
    temp = SOL_q5_T;
    SOL_q5_T = temperature(A(1), B(1), C(1), Pjsat);
    
end

SOL_q5_x = comp;

fprintf("\nQ5 Temp: " + SOL_q5_T + "\n");
fprintf("Q5 Pressure: " + SOL_q5_P + "\n");
fprintf("Q5 Liquid Comp: ");
disp(SOL_q5_x);
fprintf("Q5 Vapor Comp: ");
disp(SOL_q5_y);
    
%% %%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%

% find temp using Antoine's Equation
% returns temp in C, P in kPa
function temp = temperature(A, B, C, P)
    temp = B / (A - log(P)) - C;
end

% find saturated pressure (kPa) using Antoine's Equation temperature in Celcius
function pressure = pressure(A, B, C, T)
    pressure = exp(A - B / (T + C));
end

% Are dV, each dxi, each dyi less than ERROR
% pass the deltas into function
function bool = VLE_Condition(dV, dx, dy)
    bool = true;
    ERROR = 1e-9;
    for i = 1:size(dx)
        if (dx(i) > ERROR) || (dy(i) > ERROR)
            bool = false;
        end
    end
    
    if (dV > ERROR)
        bool = false;
    end
end

% Are each activity coeff deltas less than ERROR in consecutive iterations
function bool = gammaCondition(gamma, temp)
    bool = true;
    ERROR = 1e-9;
    for i = 1:size(gamma)
        if (abs(temp(i) - gamma(i)) > ERROR)
            bool = false;
        end
    end
end

%% phi calculations per species
% pass T in Kelvin
function phi = Phi(comp, antoine, P, T)
    
    % construct B parameter matrices

    w = [0.345, 0, 0, 0;
        0, 0.645, 0, 0;
        0, 0, 0.307, 0;
        0, 0, 0, 0.594];

    Tc = [647.1, 0, 0, 0;
        0, 513.9, 0, 0;
        0, 0, 508.2, 0;
        0, 0, 0, 563.1];

    Pc = [220.55, 0, 0, 0;
        0, 61.48, 0, 0;
        0, 0, 47.01, 0;
        0, 0, 0, 44.23];

    Zc = [0.229, 0, 0, 0;
        0, 0.240, 0, 0;
        0, 0, 0.233, 0;
        0, 0, 0, 0.260];

    Vc = [55.9, 0, 0, 0;
        0, 167.0, 0, 0;
        0, 0, 209.0, 0;
        0, 0, 0, 275.0];

    k = 0; % for essentially all species
    species = size(comp, 2);
    MAX = species;
    R = 83.145; % in bar units
    SUM = -1;
    
    P_sat = zeros(1, species);
    phi = zeros(1, species);
    delta = zeros(MAX, MAX);
    w_phi = zeros(MAX, MAX);
    Tc_phi = zeros(MAX, MAX);
    Pc_phi = zeros(MAX, MAX);
    Zc_phi = zeros(MAX, MAX);
    Vc_phi = zeros(MAX, MAX);
    B0 = zeros(MAX, MAX);
    B1 = zeros(MAX, MAX);
    B = zeros(MAX, MAX);
    B_hat = zeros(MAX, MAX);

    % Pitzer correlations and combination rule
    for i = 1:MAX
        for j = 1:MAX
            w_phi(i, j) = (w(i, i) + w(j, j)) / 2;
            Tc_phi(i, j) = ((Tc(i, i) * Tc(j, j)) ^ 0.5) * (1 - k);
            Zc_phi(i, j) = (Zc(i, i) + Zc(j, j)) / 2;
            Vc_phi(i, j) = ((Vc(i, i)^(1/3) + Vc(j, j)^(1/3))/2) ^ 3;
            Pc_phi(i, j) = Zc_phi(i, j) * R * Tc_phi(i, j) / Vc_phi(i, j);
            
            % B and delta matrices (MAX x MAX)
            B0(i, j) = 0.083 - 0.422/((T/Tc_phi(i, j)) ^ 1.6); % check for species dimensions
            B1(i, j) = 0.139 - 0.172/((T/Tc_phi(i, j)) ^ 4.2);
            B_hat(i, j) = B0(i, j) + w_phi(i, j)*B1(i, j);
            B(i, j) = B_hat(i, j) * R * Tc_phi(i, j) / Pc_phi(i, j);
        end
    end
    
    for i = 1:MAX
        for j = 1:MAX
            delta(i, j) = 2*B(i, j) - B(i, i) - B(j, j);
        end
    end

    % find saturated pressure at given T
    for i = 1:MAX
        %fprintf(antoine(i, 1) + " " + antoine(i, 2) + " " + antoine(i, 3) + " " + (T - 273.15));
        P_sat(i) = pressure(antoine(i, 1), antoine(i, 2), antoine(i, 3), T - 273.15);
    end
    
    P = P/100; % change to bar units? 
    P_sat = P_sat ./ 100;
    
    % find phi for all species
    for i = 1:species
        SUM = 0;
        for j = 1:species
            for k = 1:species
                SUM = SUM + (comp(j)*comp(k) * (2*delta(j, i) - delta(j, k)));
            end
        end
        phi(i) = exp( (B(i, i)*(P - P_sat(i)) + 0.5*P*SUM) / (R*T));
    end
    
    % detect if phi values are invalid
    
    invalid = false;
    for i = 1:MAX
        if(phi(i) > 1 || phi(i) < 0)
            invalid = true;
        end
    end
    
    %{
    if invalid
        fprintf("\nTemp: " + T);
        fprintf("\nComp: ");
        disp(comp);
        fprintf("\nPhi: ");
        disp(phi);
    end
    %}
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNI-FAC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pass T in Celcius
function activities = activity(comp, T)

    MAX = size(comp, 2);
    T = T + 273.15;
    % water(1), ethanol(2), acetone(3), 1-butanol(4)
    water = 1;
    ethanol = 2;
    acetone = 3;
    butanol = 4;

    %           Rk     Qk
    activity = [0.920, 1.400; % H2O group 7
                0.9011, 0.848; % CH3 group 1
                0.6744, 0.540; % CH2 group 1
                1.6724, 1.488; % CH3CO group 9
                1.0000, 1.200]; % OH group 5

    % subgroups
    subgroups = 5;
    H2O = [0.920, 1.400];
    CH3 = [0.0911, 0.848];
    CH2 = [0.6744, 0.540];
    CH3CO = [1.6724, 1.488];
    OH = [1.0000, 1.200];
    
    % subgroup count coefficient per species
    %    H2O CH3 CH2 CH3CO OH
    v = [1 0 0 0 0;  % H20
         0 1 1 0 1;  % CH3-CH2-OH
         0 1 0 1 0;  % CH3CO-CH3
         0 1 3 0 1]; % CH3-CH2-CH2-CH2-OH
     
    %r = [H2O(1), CH3(1)+CH2(1)+OH(1), CH3CO(1)+CH3(1), CH3(1)+3*CH2(1)+OH(1)];
    %q = [H2O(2), CH3(2)+CH2(2)+OH(2), CH3CO(2)+CH3(2), CH3(2)+3*CH2(2)+OH(2)];
    
    for i = 1:MAX
        r(i) = dot(v(i,:), activity(:,1));
        q(i) = dot(v(i,:), activity(:,2));
    end
    
    % e(k, i) (subgroup, species)
    % species     water            ethanol            acetone           butanol
    e = [H2O(2)/q(water),   0,                   0,                     0; % group H2O
         0,                 1*CH3(2)/q(ethanol), 1*CH3(2)/q(acetone),   1*CH3(2)/q(butanol); % group CH3
         0,                 1*CH2(2)/q(ethanol), 0,                     3*CH2(2)/q(butanol); % group CH2
         0,                 0,                   1*CH3CO(2)/q(acetone), 0; % CH3CO group 9
         0,                 1*OH(2)/q(ethanol),  0,                     1*OH(2)/q(butanol)]; % OH group 5

    % interaction paramters coeff matrix for subgroups
    %    H2O,     CH3,   CH2,   CH3CO,      OH
    a = [0,       300,   300, -195.40, -229.10;
         1318,      0,     0,  476.40,  986.50;
         1318,      0,     0,  476.40,  986.50;
         353.5, 156.4, 156.4,       0,  164.50;
         353.5, 156.4, 156.4,      84,       0];
     
    a = [0,       300,   300, -195.40, -229.10;
         1318,      0,     0,  476.40,  986.50;
         1318,      0,     0,  476.40,  986.50;
         472.5, 26.76, 26.76,       0,  164.50;
         353.5, 156.4, 156.4,      84,       0];

    tau = exp(-a./T);
    beta = zeros(MAX, subgroups);

    % beta 4x5 (species x subgroups)
    for i = 1:size(beta, 1)
        for k = 1:size(beta, 2)
            beta(i, k) = sum(e(:, i) .* tau(:, k));
        end
    end

    % theta 1x5 (1 x subgroups)
    theta = zeros(1, subgroups);
    for k = 1:size(theta, 2)
        theta(k) = sum(comp .* q .* e(k, :)) / sum(comp .* q);
    end

    s = zeros(1, subgroups);
    for k = 1:size(s, 2)
        s(k) = sum(theta .* tau(:, k)');
    end

    J = r ./ sum(r.*comp);
    L = q ./ sum(q.*comp);

    gammaC = 1 - J + log(J) - 5 * q.*(1 - J./L + log(J./L));
    gammaR = zeros(1, MAX);

    for i = 1:size(gammaR, 2)
        gammaR(i) = q(i) * (1 - sum(theta.*beta(i, :)./s - (e(:, i)').*log(beta(i, :)./s)));
    end

    activities = exp(gammaC + gammaR);
end