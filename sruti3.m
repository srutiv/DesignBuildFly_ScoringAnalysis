%aerodynamics: input lt, tt, vs; iterate: mt, b<5, S, power, thrust, m_prop
%output: banner length, vc, n laps
%use the outputs for scoring
Q = table;
iter = 0; %iteration count

for mt = 1:0.5:5
    for b = 0.1:0.3:1.524
        for P = 300:100:1500
            for T = 15:10:75
                for xl = 0.254:0.127:1.524 %banner length in m; minumum: 10 inches = 0.254m, max 5 feet?
                    [S, vcM2,vcM3, vt, lapsM2, lapsM3, flyM2,flyM3, M2, M3, peeps] = calculate_values(mt,b,P,T,xl);
                    if flyM2 == 0 || flyM3 == 0
                        continue
                    end
                    
                    iter = iter + 1;
                    
                    total_score = M2 + M3;
                    
                    new = table(iter,mt,b,P,T, S, vcM2, vcM3, vt,xl,peeps, M2, M3, total_score);
                    Q = [Q; new];
                    
                    %fprintf('iter = %.3f, mt = %.3f, b = %.3f, P = %.3f, T = %.3f, S = %.3f, vcM2 = %.3f, vcM3 = %.3f, vtM2 = %.3f, vtM3 = %.3f, xl = %.3f, peeps = %.3f, M2 = %.3f, M3 = %.3f \n', ...
                           %iter, mt, b, P, T, S, vc, vt, M3)
                       
                end
            end
        end
    end
end

idx = height(Q);
[~,idx] = max(Q.M3);
max_M3 = Q(idx,:)

[~,idx] = max(Q.M2);
max_M2 = Q(idx,:)

[~,idx] = max(Q.total_score);
max_total_score = Q(idx,:)


function [S, vcM2,vcM3, vt, lapsM2, lapsM3, flyM2,flyM3, M2, M3, peeps] = calculate_values(mt,b,P,T,xl)
%constants/fixed parameters
%S = area, vc = cruise velocity, x = banner length, laps = number of laps

e = 0.95; %Oswald spanwise efficiency
rho = 1.180; %density in wichita kg/m^3
g = 9.81; %gravitational acceleration in m/s^2
nu = 2; %load factor
Clmax = 1.44; %max coefficient of lift
Cd0 = 0.003; %zero lift coefficient of drag
Cd = 0.0001; %Cd0 + Clmax^2/(pi*AR*e);
lt = 6.096; %takeoff distance 20ft in m
mu = 17.97E-6; %dynamic viscosity of air; Wichita at averge 62deg F
mu_roll = 0.02; %rolling friction during taxi
f = 1.3; %factor of safety vt = fvs; otherwise, plane cannot takeoff
vmax = 14.67; % CHANGEmaximum airspeed in Wichita; used for banner Cf calculation;

%drag
Clc = 0.01; %idk how to do this because dep. on velocity
Cdc = Cd0 + (Clc^2/(pi*0.9*0.8)); %cruise Cd
Cdmax = Cd0 + (Clmax^2/(pi*40*0.8));

%prop shit
%mu_bat = 0.023;
%eff_motor = 0.4; %overall propulsive effiiency
m_bat = 0.17803 ; %battery weight in kg; Venom Lipo
%voltage = 9; %total battery voltage; Venom Lipo
%I = 42; %current in amps; Venom Lipo
% m_bat = (mu_bat*P)/(10.8*eff_motor); %battery mass
m_mot = 0.0425 ; %motor weight
m_prop = m_mot + m_bat; %total propulsion system mass


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% takeoff same for M2 and M3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = (-mt/(lt*Cd*rho))*log(1-((f^2*mt*g*Cd)/(Clmax*T))); %takeoff distance, m, t
AR = b^2/S; %aspect ratio

A = sqrt(2*T/Cd*rho*S);
B  = sqrt((T*Cd*rho*S)/(2*mt^2));
tk = (1/B)*acosh(exp(lt*(B/A))); %takeoff time

vt = A*tanh(B*tk); %takeoff velocity
vs = vt/f; %stall speed

%banner and payload
Rex = rho*vmax*xl/mu; %Reynold's number experienced by banner
Cf = 0.664/sqrt(Rex); %estimated with Blassius solution
xh = xl/5; %banner must have minimum aspect ratio of 5
Cdb = (xl*xh)*Cf/S; %banner drag; 
t = 0.0035; %thickness of banner; used 1/8" ribbon
rho_banner = 1540; %density of cotton ribbon; in kg/m^3
m_banner = xl*xh*t*rho_banner; %mass of banner

mo = 1.2*(-5.871+.8538*(S+.5976)+.03113*(S+13.12).^2); %poly fit for structural weight; input ft^2, output lb %mo = 0.454 * mo_imp; % structural weight; in kg
m_pay = (mt - m_prop - mo); %mass ALLOCATED for banner and payload
peeps = (m_pay-m_banner)/0.141748; %number of sets of passenger + luggage; 0.141748 = 5 ounces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cruise M3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cruise velocity when propulsion power is equal to drag power
%(Cd0 + Cdb)*vc^4 - (2*P/(rho*S))*vc + ((4*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%(Cd0 + Cdb)*vb^4 - (2*P/(rho*S))*vb + ((4*nu^2*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%use numerical root finder
vc_matM3 = [1.44+Cdb 0 0 ((-2*P)/(rho*S)) ((4*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here?
vb_matM3 = [1.44+Cdb 0 0 ((-2*P)/(rho*S)) ((4*nu^2*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here?

vcM3 = max(real(roots(vc_matM3)));
vbM3 = max(real(roots(vb_matM3)));

tlM3 = (2000/vcM3) + ((2*pi*mt)/(Clmax*rho*5*vbM3)); %lap time function of power, total mass
lapsM3 = ceil(600/tlM3); %actual number laps the plane can take in 10 minutes; optimal plane x=N

flyM3 = 1;
if vt < vcM3
    flyM3 = 0;
end

M3 = 2 + (lapsM3*xl)/(1.524*20); %M3 score, we want the highest value possible

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cruise M2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cruise velocity when propulsion power is equal to drag power
%(Cd0 + Cdb)*vc^4 - (2*P/(rho*S))*vc + ((4*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%(Cd0 + Cdb)*vb^4 - (2*P/(rho*S))*vb + ((4*nu^2*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%use numerical root finder
vc_matM2 = [1.44 0 0 ((-2*P)/(rho*S)) ((4*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here?
vb_matM2 = [1.44 0 0 ((-2*P)/(rho*S)) ((4*nu^2*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here?

vcM2 = max(real(roots(vc_matM2)));
vbM2 = max(real(roots(vb_matM2)));

tlM2 = (2000/vcM2) + ((2*pi*mt)/(Clmax*rho*5*vbM2)); %lap time function of power, total mass
lapsM2 = ceil(300/tlM2); %actual number laps the plane can take in 5 minutes; optimal plane x=N

flyM2 = 1;
if vt < vcM2 && lapsM2 < 3
    flyM2 = 0;
end

M2 = 1 + (peeps * vcM2)/(50*50); %FIX the denominator M3 score, we want the highest value possible

end
              
    