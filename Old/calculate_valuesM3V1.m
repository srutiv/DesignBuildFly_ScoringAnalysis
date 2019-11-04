function [S, vcM3, vtM3, lapsM3, flyM3, M3] = calculate_valuesM3(mt,b,P,T,xl)
%constants/fixed parameters
%S = area, vc = cruise velocity, x = banner length, laps = number of laps
global p foil 

opt = optimset('Display', 'off'); %Turn off warnings from fzero
e = p.e; %Oswald spanwise efficiency
rho = p.rho; %density in wichita kg/m^3
g = p.g; %gravitational acceleration in m/s^2
nu = p.nu; %load factor
Clmax = foil.Clmax; %max coefficient of lift
Cd0t = foil.Cd0t; %zero lift coefficient of drag at takeoff
Cd0c = foil.Cd0c; %zero lift coefficient of drag at cruise
Cl0c = foil.Cl0c;
lt = p.lt; %takeoff distance 20ft in m
mu = p.mu; %dynamic viscosity of air; Wichita at averge 62deg F
mu_roll = p.mu_roll; %rolling friction during taxi
f = p.f; %factor of safety vt = fvs; otherwise, plane cannot takeoff
vmax = p.vmax; % CHANGEmaximum airspeed in Wichita; used for banner Cf calculation;

%drag
Cdc = Cd0c + (Cl0c^2/(pi*0.9*0.8)); %cruise Cd; check

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
S2 = @(Cd)(pi*b^2*e*(Cd - Cd0t)/Clmax^2);
Cd = fzero(@(Cd) S1(Cd)-S2(Cd), 0.1, opt);
S = S1(Cd);
if (isnan(S))
    return;
end

% syms Cd S
% eqns = [ Cd - Cd0t - Clmax^2/(pi*(b^2/S)*e) == 0, ...
%     S - (-mt/(lt*Cd*rho))*log(1-((f^2*mt*g*Cd)/(Clmax*T))) == 0];
% K = solve(eqns, [Cd S]);
% Cd = double(K.Cd);
% S = double(K.S);

% S = (-mt/(lt*Cd*rho))*log(1-((f^2*mt*g*Cd)/(Clmax*T))); %takeoff distance, m, t; has to have same S as M2
AR = b^2/S; %aspect ratio
Cdb = (xl*xh)*Cf/S; %banner drag; 

%mo = 1.2*(-5.871+.8538*(S+.5976)+.03113*(S+13.12).^2); %poly fit for structural weight; input ft^2, output lb %mo = 0.454 * mo_imp; % structural weight; in kg
mo = 0.3*mt;
m_pay = (mt - m_prop - mo); %mass ALLOCATED for banner and payload

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
vc_matM3 = [Cdc+Cdb 0 0 ((-2*P)/(rho*S)) ((4*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here? used to be Cd0c
vb_matM3 = [Cdc+Cdb 0 0 ((-2*P)/(rho*S)) ((4*nu^2*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here? used to be Cd0c

rc = roots(vc_matM3);
rb = roots(vb_matM3);
vcM3 = max(rc(imag(rc)==0));
vbM3 = max(rb(imag(rb)==0));
% vcM3 = max(real(roots(vc_matM3)));
% vbM3 = max(real(roots(vb_matM3)));

tlM3 = (2000/vcM3) + ((2*pi*mt)/(Clmax*rho*5*vbM3)); %lap time function of power, total mass
lapsM3 = ceil(600/tlM3); %actual number laps the plane can take in 10 minutes; optimal plane x=N

flyM3 = 1;
% if vtM3*3 < vcM3
%     flyM3 = 0;
% end

if vcM3 < 1.3*vsM3 || m_pay < 0
    flyM3 = 0;
end

M3 = 2 + (lapsM3*xl)/(1.524*20); %M3 score, we want the highest value possible

end