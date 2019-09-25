%aerodynamics: input lt, tt, vs; iterate: mt, b<5, S, power, thrust, m_prop
%output: banner length, vc, n laps
%use the outputs for scoring

for mt = 5:5
    for b = 0.20
        for P = 500
            for T = 20
                [S, vc, x, laps] = calculate_values(mt,b,P,T);
            end
        end
    end
end

function [S, vc, x, laps] = calculate_values(mt,b,P,T)
%constants/fixed parameters
%S = area, vc = cruise velocity, x = banner length, laps = number of laps
eff_motor = 0.4;
e = 0.1; %efficiency factor; i picked a rando number
rho = 1.180; %density in wichita kg/m^3
Clmax = 1.4;
Clc = 0.8; %cruise Cl
Cd0 = 0.01;
vs = 10; %stall speed in m/s
%tt = 10*60; %total flight time in seconds
lt = 6.096; %takeoff distance 20ft in m
g = 9.81; %gravitational acceleration in m/s^2
nu = 2; %load factor
x = 1;
%where the fuck does T factor in?
%where does lt factor in?

S = (2*mt*g)/(Clmax*rho*vs^2); %function of mt and vs; in m/s

mu_bat = 0.023;
m_bat = (mu_bat*P)/(10.8*eff_motor); %battery mass
mprop = 0.5 + m_bat; %total propulsion system mass
mo = 1.2*(-5.871+.8538*(S+.5976)+.03113*(S+13.12).^2); %poly fit for structural weight; input ft^2, output lb %mo = 0.454 * mo_imp; % structural weight; in kg
mpay = (mt - mprop - mo); %mass ALLOCATED for banner and payload

Cdb = 0.0001; %banner drag; rando
Cdc = Cd0 + (Clc^2/(pi*0.9*0.8)) + x*Cdb; %cruise Cd; Cdb = drag of banner per unit length
Cdmax = Cd0 + (Clmax^2/(pi*40*0.8)); % + (N*(Cdo/S));
%incorporate banner somehow into drag (which will drive total drag and vc)

%cruise velocity when propulsion power is equal to drag power
%Cd0*vc^4 - (2*P/(rho*S))*vc + ((4*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%Cd0*vb^4 - (2*P/(rho*S))*vb + ((4*nu^2*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%use numerical root finder
vc_mat = [Cd0 0 0 ((-2*P)/(rho*S)) ((4*mt^2*g^2)/(rho^2*S*b^2*pi*e))];
vb_mat = [Cd0 0 0 ((-2*P)/(rho*S)) ((4*nu^2*mt^2*g^2)/(rho^2*S*b^2*pi*e))];

vc = roots(vc_mat);
vb = roots(vb_mat);

vc = vc(1);
vb = vb(1);

tt = ((2000)/vc) + ((2*pi*mt)/(Clmax*rho*5*vb)); %total lap time function of power, total mass
laps = ceil(600/tt); %actual number laps the plane can take in 10 minutes; optimal plane x=N

end
                
                
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
% %everything a function of total mass, takeoff speed, and power
% %find the fastest plane that can carry the most payloads
% %#laps should equal #payloads
% %THIS IS SUPER ROUGH; NEED TO: factor in payload drag, recheck structural
% %weight poly fit
% 
% clc
% clear all
% 
% %iteration ranges
% m_range = 1:1:13; %range of total mass in kg
% P_range = 100:1000:5100; %range of power in watts
% v_range = 1:1:12; %range of takeoff velocity in m/s
% score = {['m,  ', 'P,  ', 'v,  ', 'S,  ', 'N,  ', 'x']}; 
% 
% for m = m_range
%     for P = P_range
%         for v = v_range
%             [mn, mp, S, N, x] = calculate_values(m,P,v);
%             i = length(score) + 1;
%             if N > 0 && N == x %plane must carry at least 1 payload and number of payloads == number of laps
%                 score{i} = [m, P, v, S, N, x];
%             end
%         end
%     end
% end
% 
% score = score';
% 
% function [mpay, mprop, S, N, x] = calculate_values(m,P,v)
% %constants/fixed parameters
% eff_motor = 0.4;
% rho = 1.180; %density in wichita kg/m^3
% Clmax = 1.4;
% Clc = 0.8; %cruise Cl
% Cdo = 0.01;
% 
% S = (2.4*m*9.81)/(Clmax*rho*v^2); %function of total mass and takeoff velocity; output in m^2 %S_imp = 10.764 * S; %function of mt and v, output in ft^2
% mu_bat = 0.023;
% m_bat = (mu_bat*P)/(10.8*eff_motor); %battery mass
% mprop = 0.5 + m_bat; %total propulsion system mass
% mo = 1.2*(-5.871+.8538*(S+.5976)+.03113*(S+13.12).^2); %poly fit for structural weight; input ft^2, output lb %mo = 0.454 * mo_imp; % structural weight; in kg
% mpay = (m - mprop - mo); %mass ALLOCATED for payload
% N = ceil(mpay/0.086); %number of payloads
% 
% Cdc = Cdo + (Clc^2/(pi*0.9*0.8)) + (N*(Cdo/S)); %cruise Cd
% vc = ((Clmax*P*v^2)/(1.2*Cdc*m*9.81))^(1/3);
% Cdmax = Cdo + (Clmax^2/(pi*40*0.8)) + (N*(Cdo/S));
% 
% vb = ((Clmax*P*v^2)/(1.2*Cdmax*m*9.81))^(1/3); %turn velocity
% tl = ((2000)/vc) + ((2*pi*m)/(Clmax*rho*5*vb)); %lap time function of power, total mass
% x = ceil(600/tl); %actual number laps the plane can take in 10 minutes; optimal plane x=N
% end
% 
% nu = 2; %max load factor; n= L/W
% L = 2*9.81*mt; %max lift generated; during turn radius
                
    