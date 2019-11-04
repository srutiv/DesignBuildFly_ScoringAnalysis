%aerodynamics: input lt, tt, vs; iterate: mt, b<5, S, power, thrust, m_prop
%output: banner length, vc, n laps
%use the outputs for scoring
tic
warning off;
Q = table;
%iter = 0; %iteration count
global p foil
p = getConstants();
foil = get_Airfoil('mh32_200000.txt', 'mh32_500000.txt');
fprintf('loaded airfoil data \n');

iter = 0;
for mt = 1:0.5:8 %don't really want to build a heavy-ass plane
    for b = 0.5:0.1:1.524 %maximum wingspan = 5ft
        for P = 300:100:600 %range of values provided by Will; range limited by budget, safety, and team's comfort and experience
            for T = 25:5:50 %range of values provided by Will
                for xl = 0.254:0.1:1.524 %banner length in m; minumum: 10 inches = 0.254m, max 5 feet?
                    try
                        [SM2, vcM2, vtM2, lapsM2, flyM2, peeps, prod2,tt2, profileM2, mStruct] = calculate_valuesM2(mt,b,P,T,xl);
                        [SM3, vcM3, vtM3, lapsM3, flyM3, prod3, profileM3] = calculate_valuesM3(mt,b,P,T,xl, mStruct);
                    catch
                        continue;
                    end
                    
                    if flyM2 == 0 || flyM3 == 0
                        continue
                    end
                    
                    iter = iter + 1;
                    
                    %new = table(iter,mt,b,P,T, S, vcM2, vcM3, vt,xl,peeps, M2, M3, total_score);
                    new = table(iter,mt,b,P,T,SM2,SM3,vcM2,vtM2,lapsM2,xl,peeps,vcM3,vtM3,lapsM3,tt2, prod2, prod3, profileM2, profileM3);
                    Q = [Q; new];
                    
                       
                end
            end
        end
    end
end

idx = height(Q);
[~,idx] = max(Q.prod2);
max_M2 = Q(idx,:);
max_M2 = max_M2(1,'prod2');
max_M2 = table2array(max_M2);

[~,idx] = max(Q.prod3);
max_M3 = Q(idx,:);
max_M3 = max_M3(1,'prod3');
max_M3 = table2array(max_M3);

A = table2cell(Q);
[m,n] = size(A);
NewCol = rand(m,1);
%A = [A NewCol];

for i = 1:m
    num = A{i,12}; denom = A{i,16};
    M2 = 1 + (num/denom)/max_M2; %(peeps/tt2)/max_M2;
    num = A{i,15}; denom = A{i,11};
    M3 = 2 + (num*denom)/max_M3; %(lapsM3*xl)/max_M3;
    total_score = M2 + M3;
    A{i,n+1} = total_score;
end

W = cell2table(A,'VariableNames',{'iter','mt','b','P','T','SM2', ...
    'SM3','vcM2','vtM2','lapsM2','xl','peeps','vcM3','vtM3','lapsM3','tt2', 'prod2', 'prod3','profileM2', 'profileM3', 'total_score',});

% [~,idx] = max(W.prod2);
% max_prod2 = W(idx,:);
% disp(max_prod2.profileM2)
% disp(max_prod2.profileM3)
% 
% [~,idx] = max(W.prod3);
% max_prod3 = W(idx,:);
% disp(max_prod3.profileM2)
% disp(max_prod3.profileM3)

[~,idx] = max(W.total_score);
max_total_score = W(idx,:);

fprintf('plane that will give the highest total score is:   \n')
fprintf('one with a M2 profile of:    \n')
disp(max_total_score.profileM2)

fprintf('and one with a M3 profile of:    \n')
disp(max_total_score.profileM3)


ordered = sortrows(W,21,'descend');

toc


function p = getConstants()
    global p
    p.e = 0.95; %Oswald spanwise efficiency
    p.rho = 1.180; %density in wichita kg/m^3
    p.g = 9.81; %gravitational acceleration in m/s^2
    p.nu = 2; %load factor
    p.lt = 6.096; %takeoff distance 20ft in m
    p.mu = 17.97E-6; %dynamic viscosity of air; Wichita at averge 62deg F
    p.mu_roll = 0.2; %rolling friction during taxi
    p.f = 1.3; %factor of safety vt = fvs; otherwise, plane cannot takeoff
    p.vmax = 14.67; % CHANGEmaximum airspeed in Wichita; used for banner Cf calculation;
    p.mu_bat = (0.490/4); %mass/cell for 4s battery %change to make dependent on how many no cells
    p.eta = 0.6; %mechanical efficiency factor
    p.nom_volt = 3.7; %in volts; nominal voltage for lipos
    p.capacity = 5000; %battery capacity in mAh
    p.I_pack = (p.capacity*10^-3)/(5/60); %current draw of pack in Amps; 10/60 = mission time in hours
    %m_bat = 0.17803 ; %battery weight in kg; Venom Lipo
    %voltage = 9; %total battery voltage; Venom Lipo
    %I = 42; %current in amps; Venom Lipo
    p.m_mot = 0.3 ; %upper limit for motor weight in kg; fixed

end

function foil = get_Airfoil(airofil_takeoff, airfoil_cruise)
    global foil
    
    % Takeoff, Re = 200000
    file = textread(airofil_takeoff, '%s', 'delimiter', '\n','whitespace', ' ');
    i = 12; % Find where the data begins

    % Store all results in arrays
    alphas = []; CLs = []; CDs = []; CDps = [];
    if(i ~= length(file) && ~isempty(file(i+1)))
        j = i+1;
        while (j<=length(file) && ~isempty(file(j)))
            results = textscan(char(file(j)), '%f');
            alphas(j-i) = results{1}(1);
            CLs(j-i) = results{1}(2);
            CDs(j-i) = results{1}(3);
            CDps(j-i) = results{1}(4);
            j = j+1;
        end
    end

    Cl0t=min(abs(CLs)); % smallest cl in file (closest to zero lift)
    foil.Cd0t=CDs(abs(CLs)==Cl0t); % Cd at zero lift; used for vt calculations
    foil.Clmax = max(CLs); %2d wing maximum CL
    
    % Cruise, Re = 500000;
    file = textread(airfoil_cruise, '%s', 'delimiter', '\n','whitespace', ' ');
    i = 12; % Find where the data begins

    % Store all results in arrays
    alphas = []; CLs = []; CDs = []; CDps = [];
    if(i ~= length(file) && ~isempty(file(i+1)))
        j = i+1;
        while (j<=length(file) && ~isempty(file(j)))
            results = textscan(char(file(j)), '%f');
            alphas(j-i) = results{1}(1);
            CLs(j-i) = results{1}(2);
            CDs(j-i) = results{1}(3);
            CDps(j-i) = results{1}(4);
            j = j+1;
        end
    end

    Cl0c=min(abs(CLs)); % smallest cl in file (closest to zero lift)
    foil.Cl0c = Cl0c;
    foil.Cd0c = CDs(abs(CLs)==Cl0c); %zero lift coefficient of drag at cruise; used to solve for vc and vb

end