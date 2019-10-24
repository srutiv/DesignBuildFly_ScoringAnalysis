%aerodynamics: input lt, tt, vs; iterate: mt, b<5, S, power, thrust, m_prop
%output: banner length, vc, n laps
%use the outputs for scoring
warning off;
Q = table;
iter = 0; %iteration count
p = getConstants();
foil = get_Airfoil('mh32_200000.txt', 'mh32_500000.txt');

for mt = 1:0.5:3
    for b = 0.5:0.1:1.524
        for P = 300:20:600
            for T = 25 %20:5:50
                for xl = 0.254:0.127:1 %1.524 %banner length in m; minumum: 10 inches = 0.254m, max 5 feet?
%                     fprintf('%f %f %f %f %f\n', mt, b, P, T, xl);
                    SM2 = 0; vcM2 = 0; vtM2 = 0; lapsM2 = 0; flyM2 = 0; flyM3 = 0; M2 = 0; M3 = 0; peeps = 0;
                    try
                        [SM2, vcM2, vtM2, lapsM2, flyM2, M2, peeps] = calculate_valuesM2(mt,b,P,T,xl,p);
                        [SM3, vcM3, vtM3, lapsM3, flyM3, M3] = calculate_valuesM3(mt,b,P,T,xl,p);
                    catch
                        continue;
                    end
                    if flyM2 == 0 || flyM3 == 0
                        continue
                    end
                    iter = iter + 1;
                    
                    total_score = M2 + M3;
                    
                    %new = table(iter,mt,b,P,T, S, vcM2, vcM3, vt,xl,peeps, M2, M3, total_score);
                    new = table(iter,mt,b,P,T,SM2,SM3, vcM2,vtM2,lapsM2,xl,peeps,vcM3,vtM3,lapsM3,M2,M3,total_score);
                    Q = [Q; new];
                    
                    %fprintf('iter = %.3f, mt = %.3f, b = %.3f, P = %.3f, T = %.3f, S = %.3f, vcM2 = %.3f, vcM3 = %.3f, vtM2 = %.3f, vtM3 = %.3f, xl = %.3f, peeps = %.3f, M2 = %.3f, M3 = %.3f \n', ...
                           %iter, mt, b, P, T, S, vc, vt, M3)
                       
                end
            end
        end
    end
end

idx = height(Q);
[~,idx] = max(Q.M2);
max_M2 = Q(idx,:)

[~,idx] = max(Q.M3);
max_M3 = Q(idx,:)

[~,idx] = max(Q.total_score);
max_total_score = Q(idx,:)


function [S, vcM2,vtM2,lapsM2,flyM2,M2,peeps] = calculate_valuesM2(mt,b,P,T,xl,p)
%constants/fixed parameters
%S = area, vc = cruise velocity, x = banner length, laps = number of laps

opt = optimset('Display', 'off'); %Turn of warnings from fzero
e = p.e; %Oswald spanwise efficiency
rho = p.rho; %density in wichita kg/m^3
g = p.g; %gravitational acceleration in m/s^2
nu = p.nu; %load factor
Clmax = p.Clmax; %max coefficient of lift
Cd0 = p.Cd0; %zero lift coefficient of drag
lt = p.lt; %takeoff distance 20ft in m
mu = p.mu; %dynamic viscosity of air; Wichita at averge 62deg F
mu_roll = p.mu_roll; %rolling friction during taxi
f = p.f; %factor of safety vt = fvs; otherwise, plane cannot takeoff
vmax = p.vmax; % CHANGEmaximum airspeed in Wichita; used for banner Cf calculation;

%drag
Clc = 0.01; %idk how to do this because dep. on velocity
Cdc = Cd0 + (Clc^2/(pi*0.9*0.8)); %cruise Cd
Cdmax = Cd0 + (Clmax^2/(pi*40*0.8));

%prop
mu_bat = p.mu_bat; %mass/cell
eta = p.eta; %mechanical efficiency factor
nom_volt = p.nom_volt; %in volts; nominal voltage for lipos
capacity = p.capacity; %battery capacity in mAh
I_pack = p.I_pack; %current draw of pack in Amps; 10/60 = mission time in hours
%m_bat = 0.17803 ; %battery weight in kg; Venom Lipo
%voltage = 9; %total battery voltage; Venom Lipo
%I = 42; %current in amps; Venom Lipo
m_mot = p.m_mot; %upper limit for motor weight in kg; fixed
m_prop = m_mot + ((mu_bat*P)/(nom_volt*I_pack*eta)); %total propulsion system mass


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% takeoff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S1 = @(Cd)(-mt./(lt*Cd*rho)).*log(1-((f^2*mt*g*Cd)/(Clmax*T)));
S2 = @(Cd)(pi*b^2*e*(Cd - Cd0)/Clmax^2);
Cd = fzero(@(Cd) S1(Cd)-S2(Cd), 0.1, opt);
S = S1(Cd);
if (isnan(S))
    return;
end

A = sqrt(2*T/(Cd*rho*S));
B  = sqrt((T*Cd*rho*S)/(2*mt^2));
tk = (1/B)*acosh(exp(lt*(B/A))); %takeoff time

vtM2 = A*tanh(B*tk); %takeoff velocity
vsM2 = vtM2/f; %stall speed

%banner and payload
Rex = rho*vmax*xl/mu; %Reynold's number experienced by banner
Cf = 0.664/sqrt(Rex); %estimated with Blassius solution
xh = xl/5; %banner must have minimum aspect ratio of 5
Cdb = 0; %banner drag;  not engaged in M2
t = 0.0035; %thickness of banner; used 1/8" ribbon
rho_banner = 1540; %density of cotton ribbon; in kg/m^3
m_banner = xl*xh*t*rho_banner; %mass of banner

%mo = 0.453592*(1.2*(-5.871+.8538*((S*10.764)+.5976)+.03113*(S+13.12).^2)); %poly fit for structural weight; input m^2, output kg; do we trust this
mo = 0.3*mt; 
m_pay = (mt - m_prop - mo); %mass ALLOCATED for banner and payload
peeps = (m_pay-m_banner)/0.141748; %number of sets of passenger + luggage; 0.141748 = 5 ounces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cruise M2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cruise velocity when propulsion power is equal to drag power
%(Cd0 + Cdb)*vc^4 - (2*P/(rho*S))*vc + ((4*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%(Cd0 + Cdb)*vb^4 - (2*P/(rho*S))*vb + ((4*nu^2*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%use numerical root finder
Cdb = 0.561 * (xl/xh).^-0.480;
vc_matM2 = [Cd0 0 0 ((-2*P)/(rho*S)) ((4*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here?
vb_matM2 = [Cd0 0 0 ((-2*P)/(rho*S)) ((4*nu^2*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here?
rc = roots(vc_matM2);
rb = roots(vb_matM2);
vcM2 = max(rc(imag(rc)==0));
vbM2 = max(rb(imag(rb)==0));


tlM2 = (2000/vcM2) + ((2*pi*mt)/(Clmax*rho*5*vbM2)); %lap time function of power, total mass
lapsM2 = ceil(300/tlM2); %actual number laps the plane can take in 5 minutes; optimal plane x=N

flyM2 = 1;
% fprintf('%.2f vtm2\n', vtM2);
% fprintf('%.2f vcm2\n', vcM2);
% fprintf('%.2f lapsM2\n', lapsM2);
if vtM2*3 < vcM2 || lapsM2 < 3
    flyM2 = 0;
end

M2 = 1 + (peeps * vcM2)/100; %FIX the denominator M3 score, we want the highest value possible

end


function [S, vcM3, vtM3, lapsM3, flyM3, M3] = calculate_valuesM3(mt,b,P,T,xl,p)
%constants/fixed parameters
%S = area, vc = cruise velocity, x = banner length, laps = number of laps

opt = optimset('Display', 'off'); %Turn off warnings from fzero
e = p.e; %Oswald spanwise efficiency
rho = p.rho; %density in wichita kg/m^3
g = p.g; %gravitational acceleration in m/s^2
nu = p.nu; %load factor
Clmax = foil.Clmax; %max coefficient of lift
Cd0t = foil.Cd0t; %zero lift coefficient of drag at takeoff
Cd0c = foil.Cd0c; %zero lift coefficient of drag at cruise
lt = p.lt; %takeoff distance 20ft in m
mu = p.mu; %dynamic viscosity of air; Wichita at averge 62deg F
mu_roll = p.mu_roll; %rolling friction during taxi
f = p.f; %factor of safety vt = fvs; otherwise, plane cannot takeoff
vmax = p.vmax; % CHANGEmaximum airspeed in Wichita; used for banner Cf calculation;

%drag
Clc = 0.01; %idk how to do this because dep. on velocity
Cdc = Cd0c + (Clc^2/(pi*0.9*0.8)); %cruise Cd
Cdmax = Cd0c + (Clmax^2/(pi*40*0.8)); %max Cd at cruise

%prop
mu_bat = p.mu_bat; %mass/cell
eta = p.eta; %mechanical efficiency factor
nom_volt = p.nom_volt; %in volts; nominal voltage for lipos
capacity = p.capacity; %battery capacity in mAh
I_pack = p.I_pack; %current draw of pack in Amps; 10/60 = mission time in hours
%m_bat = 0.17803 ; %battery weight in kg; Venom Lipo
%voltage = 9; %total battery voltage; Venom Lipo
%I = 42; %current in amps; Venom Lipo
m_mot = p.m_mot; %upper limit for motor weight in kg; fixed
m_prop = m_mot + ((mu_bat*P)/(nom_volt*I_pack*eta)); %total propulsion system mass

%banner and payload
Rex = rho*vmax*xl/mu; %Reynold's number experienced by banner
Cf = 0.664/sqrt(Rex); %estimated with Blassius solution
xh = xl/5; %banner must have minimum aspect ratio of 5
t = 0.0035; %thickness of banner; used 1/8" ribbon
rho_banner = 1540; %density of cotton ribbon; in kg/m^3
m_banner = xl*xh*t*rho_banner; %mass of banner

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% takeoff
S1 = @(Cd)(-mt./(lt*Cd*rho)).*log(1-((f^2*mt*g*Cd)/(Clmax*T)));
S2 = @(Cd)(pi*b^2*e*(Cd - Cd0)/Clmax^2);
Cd = fzero(@(Cd) S1(Cd)-S2(Cd), 0.1, opt);
S = S1(Cd);
if (isnan(S))
    return;
end
% S = (-mt/(lt*Cd*rho))*log(1-((f^2*mt*g*Cd)/(Clmax*T))); %takeoff distance, m, t; has to have same S as M2
AR = b^2/S; %aspect ratio
Cdb = (xl*xh)*Cf/S; %banner drag; 

%mo = 1.2*(-5.871+.8538*(S+.5976)+.03113*(S+13.12).^2); %poly fit for structural weight; input ft^2, output lb %mo = 0.454 * mo_imp; % structural weight; in kg
mo = 0.3*mt;
m_pay = (mt - m_prop - mo); %mass ALLOCATED for banner and payload
peeps = (m_pay-m_banner)/0.141748; %number of sets of passenger + luggage; 0.141748 = 5 ounces

mt = mt - m_pay + m_banner; %redefine mt because no passengers for M3

A = sqrt(2*T/Cd*rho*S);
B  = sqrt((T*Cd*rho*S)/(2*mt^2));
tk = (1/B)*acosh(exp(lt*(B/A))); %takeoff time

vtM3 = A*tanh(B*tk); %takeoff velocity
vsM3 = vtM3/f; %stall speed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cruise M3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cruise velocity when propulsion power is equal to drag power
%(Cd0 + Cdb)*vc^4 - (2*P/(rho*S))*vc + ((4*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%(Cd0 + Cdb)*vb^4 - (2*P/(rho*S))*vb + ((4*nu^2*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%use numerical root finder
vc_matM3 = [Cd0c+Cdb 0 0 ((-2*P)/(rho*S)) ((4*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here?
vb_matM3 = [Cd0c+Cdb 0 0 ((-2*P)/(rho*S)) ((4*nu^2*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here?

rc = roots(vc_matM3);
rb = roots(vb_matM3);
vcM3 = max(rc(imag(rc)==0));
vbM3 = max(rb(imag(rb)==0));
% vcM3 = max(real(roots(vc_matM3)));
% vbM3 = max(real(roots(vb_matM3)));

tlM3 = (2000/vcM3) + ((2*pi*mt)/(Clmax*rho*5*vbM3)); %lap time function of power, total mass
lapsM3 = ceil(600/tlM3); %actual number laps the plane can take in 10 minutes; optimal plane x=N

flyM3 = 1;
if vtM3*3 < vcM3
    flyM3 = 0;
end

M3 = 2 + (lapsM3*xl)/(1.524*20); %M3 score, we want the highest value possible

end

function p = getConstants()
    
    p.e = 0.95; %Oswald spanwise efficiency
    p.rho = 1.180; %density in wichita kg/m^3
    p.g = 9.81; %gravitational acceleration in m/s^2
    p.nu = 2; %load factor
    p.lt = 6.096; %takeoff distance 20ft in m
    p.mu = 17.97E-6; %dynamic viscosity of air; Wichita at averge 62deg F
    p.mu_roll = 0.02; %rolling friction during taxi
    p.f = 1.3; %factor of safety vt = fvs; otherwise, plane cannot takeoff
    p.vmax = 14.67; % CHANGEmaximum airspeed in Wichita; used for banner Cf calculation;

    p.mu_bat = (0.490/4); %mass/cell for 4s battery %change to make dependent on how many no cells
    p.eta = 0.4; %mechanical efficiency factor
    p.nom_volt = 3.7; %in volts; nominal voltage for lipos
    p.capacity = 5000; %battery capacity in mAh
    p.I_pack = (p.capacity*10^-3)/(5/60); %current draw of pack in Amps; 10/60 = mission time in hours
    %m_bat = 0.17803 ; %battery weight in kg; Venom Lipo
    %voltage = 9; %total battery voltage; Venom Lipo
    %I = 42; %current in amps; Venom Lipo
    p.m_mot = 0.3 ; %upper limit for motor weight in kg; fixed

end

function foil = get_Airfoil(airofil_takeoff, airfoil_cruise)
    
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

    cl0t=min(abs(clt)); % smallest cl in file (closest to zero lift)
    foil.Cd0t=cdt(abs(clt)==cl0t); %zero lift coefficient of drag at takeoff; used to solve for vt
    foil.Clmaxt = 0.9*max(clt); %2d wing maximum; max coefficient of lift at takeoff
    
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

    cl0c=min(abs(cl)); % smallest cl in file (closest to zero lift)
    foil.Cd0c = cd(abs(cl)==cl0c); %zero lift coefficient of drag at cruise; used to solve for vc and vb

end