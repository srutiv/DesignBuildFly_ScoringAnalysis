function [SM2, vcM2, vtM2, lapsM2, flyM2, peeps, prod2,tt, profileM2, mo] = calculate_valuesM2(mt,b,P,T,xl, S_sens);
%constants/fixed parameters
global p foil

opt = optimset('Display', 'off'); %Turn of warnings from fzero
e = p.e; %Oswald spanwise efficiency
rho = p.rho; %density in wichita kg/m^3
g = p.g; %gravitational acceleration in m/s^2
nu = p.nu; %load factor
lt = p.lt; %takeoff distance 20ft in m
mu = p.mu; %dynamic viscosity of air; Wichita at averge 62deg F
mu_roll = p.mu_roll; %rolling friction during taxi
f = p.f; %factor of safety vt = fvs; otherwise, plane cannot takeoff
vmax = p.vmax; % CHANGE  maximum airspeed in Wichita; used for banner Cf calculation;
fos = p.fos; %factor of safety for area and thrust

%airfoil
Cd0c = foil.Cd0c; %zero lift coefficient of drag
Cd0t = foil.Cd0t;
Cl0c = foil.Cl0c;
Clmax = foil.Clmax; %max coefficient of lift

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
if S_sens == 0 %regular run
    S1 = @(Cd)(-mt./(lt*Cd*rho)).*log(1-((f^2*mt*g*Cd)/(Clmax*T))); %factor of safety for area and thrust)));
    S2 = @(Cd)(pi*b^2*e*(Cd - Cd0t)/Clmax^2);
    Cd = fzero(@(Cd) S1(Cd)-S2(Cd), 0.1, opt);
    S = S1(Cd);
    if (isnan(S))
        return;
    end
else %sensitivity run; S_sense nonzero 0
    %%%is this correct? did I recalculate Cd correctly?? %%%%
    
    S1 = @(Cd)(-mt./(lt*Cd*rho)).*log(1-((f^2*mt*g*Cd)/(Clmax*T)));
    S2 = @(Cd)(pi*b^2*e*(Cd - Cd0t)/Clmax^2);
    Cd = fzero(@(Cd) S1(Cd)-S2(Cd), 0.1, opt);
    S = S_sens;
    if (isnan(S))
        return;
    end
end


% syms Cd S
% eqns = [ Cd - Cd0t - Clmax^2/(pi*(b^2/S)*e) == 0, ...
%     S - (-mt/(lt*Cd*rho))*log(1-((f^2*mt*g*Cd)/(Clmax*T))) == 0];
% K = solve(eqns, [Cd S]);
% Cd = double(K.Cd);
% S = double(K.S);

%S = (-mt/(lt*Cd*rho))*log(1-((f^2*mt*g*Cd)/(Clmax*T))); %takeoff distance, m, t

AR = b^2/S; %aspect ratio

A = sqrt(2*T/(Cd*rho*S));
B  = sqrt((T*Cd*rho*S)/(2*mt^2));
tk = (1/B)*acosh(exp(lt*(B/A))); %takeoff time

vtM2 = A*tanh(B*tk); %takeoff velocity
vsM2 = vtM2/f; %stall speed

%banner
% Rex = rho*vmax*xl/mu; %Reynold's number experienced by banner
% Cf = 0.664/sqrt(Rex); %estimated with Blassius solution
% xh = xl/5; %banner must have minimum aspect ratio of 5
% Cdb = 0; %banner drag;  not engaged in M2
% t = 0.0035; %thickness of banner; used 1/8" ribbon
% rho_banner = 1540; %density of cotton ribbon; in kg/m^3
% m_banner = xl*xh*t*rho_banner; %mass of banner

%payload
%mo = 0.453592*(1.2*(-5.871+.8538*((S*10.764)+.5976)+.03113*(S+13.12).^2)); %poly fit for structural weight; input m^2, output kg; do we trust this
c = S/b;
mo = getStructuralWeight(114,1050,b*1000,c*1000); %estimation of structural weight; need a better way to do this
m_pay = (mt - m_prop - mo); %mass ALLOCATED for banner and payload
peeps = (m_pay)/0.141748; %number of sets of passenger + luggage; 0.141748 = 5 ounces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cruise M2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cruise velocity when propulsion power is equal to drag power
%(Cd0 + Cdb)*vc^4 - (2*P/(rho*S))*vc + ((4*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%(Cd0 + Cdb)*vb^4 - (2*P/(rho*S))*vb + ((4*nu^2*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%use numerical root finder
%Cdb = 0.561 * (xl/xh).^-0.480; %no banner in M2

% vc_matM2 = [Cdc 0 0 ((-2*P)/(rho*S)) ((4*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; 
% vb_matM2 = [Cdc 0 0 ((-2*P)/(rho*S)) ((4*nu^2*mt^2*g^2)/(rho^2*S*b^2*pi*e))];
% rc = roots(vc_matM2);
% rb = roots(vb_matM2);
% vcM2 = max(rc(imag(rc)==0));
% vbM2 = max(rb(imag(rb)==0));

% syms v4
% %eqns = T*vcM2 - P == 0, 
% eqns = [(0.5*rho*(v4^3)*S)*(Cd0c + (((2*mt*9.8/(S*rho*v4^2))^2)/(pi*e*(b^2/S)))) - P*eff == 0];
% J = solve(eqns, v4);
% v4 = double(max(abs(J)));


fun1 = @(vcM2) (0.5*rho*(vcM2^3)*S)*(Cd0c + (((2*mt*9.8/(S*rho*vcM2^2))^2)/(pi*e*(b^2/S)))) - P*eta; % function
x0 = 10; % initial point; possibly change
vcM2 = fzero(fun1,x0);
if (isnan(vcM2))
    return;
end

% syms vcM2
% %eqns = T*vcM2 - P == 0, 
% eqns = [(0.5*rho*(vcM2^3)*S)*(Cd0c + (((2*mt*9.8/(S*rho*vcM2^2))^2)/(pi*e*(b^2/S)))) - P*eta == 0];
% J = solve(eqns, vcM2);)
% vcM2 = double(max(abs(J)));

fun2 = @(vbM2) (0.5*rho*(vbM2^3)*S)*(Cd0c + (((2*mt*9.8*nu/(S*rho*vbM2^2))^2)/(pi*e*(b^2/S)))) - P*eta; % function
x0 = 10; % initial point; possibly change
vbM2 = fzero(fun2,x0);
if (isnan(vbM2))
    return;
end

% syms vbM2 
% eqns = [(0.5*rho*(vbM2^3)*S)*(Cd0c + (((2*mt*9.8*nu/(S*rho*vbM2^2))^2)/(pi*e*(b^2/S)))) - P*eta == 0];
% K = solve(eqns, vbM2);
% vbM2 = double(max(abs(K)));

tlM2 = (2000/vcM2) + ((2*pi*mt)/(Clmax*rho*5*vbM2)); %lap time function of power, total mass
lapsM2 = ceil(300/tlM2); %actual number laps the plane can take in 5 minutes; optimal plane x=N

SM2 = S;

flyM2 = 1;
% if vtM2*3 < vcM2 || lapsM2 < 3
%     flyM2 = 0;
% end

if vcM2 < 1.3*vsM2 || m_pay < 0
    flyM2 = 0;
end

tt = tlM2 + tk;
prod2 = peeps/tt; 

%package mission profile;
profileM2.mt2 = mt;
profileM2.b2 = b;
profileM2.P2 = P;
profileM2.T2 = T;
profileM2.xl2 = 0;

profileM2.S2 = S;
profileM2.AR2 = AR;

profileM2.m_prop2 = m_prop;
profileM2.mo2 = mo;
profileM2.m_pay2 = m_pay;

profileM2.vtM2 = vtM2;
profileM2.vsM2 = vsM2;
profileM2.vcM2 = vcM2;
profileM2.vbM2 = vbM2;


profileM2.peeps2 = peeps;
profileM2.lapsM2 = lapsM2;
profileM2.tt2 = tt;
profileM2.flyM2 = flyM2;

end