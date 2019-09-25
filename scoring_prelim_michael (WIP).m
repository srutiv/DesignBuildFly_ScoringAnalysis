% Michael Zakoworotny
% DBF
% Preliminary Scoring Testing

%Conversions
ftTOm = 0.3048; % ft to m conversion
ozTOkg = 0.0283495; %oz to kg covnersion

%First approximations of geometric values
% b = 5 * ftTOm; %m, span
t_3 = 10 * 60; %seconds, M3 duration
d_takeoff = 20 * ftTOm; %m, takeoff distance
c_root = 1 * ftTOm; %root chord
lam = 0.9; %taper ratio
e = 0.8; %Oswald efficiency, for "other aircraft" according to Raymer
CD0 = 0.02; %Zero-lift coefficient of drag, for a "clean propeller aircraft" according to Raymer
al = 0; %Angle of attack
m_pass = 4 * ozTOkg; %mass of passenger
m_lugg = 1 * ozTOkg; %mass of luggage
m_struct = 1;
m_wing = 0.5;
Cf_ban = 7.72e-3; %Cf = Cd for objects with no pressure drag
ae = aircraft(b, c_root, lam, e, CD0, al); %Store these parameters in a struct

%Lithium Ion Data, https://www.claytonpower.com/products/lithium-ion-batteries/specifications/
q = 100 * 3600; %Amp*hours -> Amp*second
V = 12.6; % Volts
I = 63; %Continuous discharge current
eta = 0.7; %efficiency of motor
P_motor = I*V; %Power of motor
% batt_dens = 10 * P_motor/1000*3600; %kg/kWh
m_bat = 178.03/1000; %kg mass of battery, very rought estimate
m_motor = 3.10 * ozTOkg; %cobra motor mass, https://www.cobramotorsusa.com/motors-2221-16.html
P_mech = eta*P_motor;

%Constants
c.rho = 1.14; %kg/m^3, calculated density in Wichita using altitude
c.visc = 1.765e-5;
c.mu = 0.02;
c.g = 9.81;

% n_pass = 1:1:20;
% v_c = 1:0.5:100;
configs = {};
i = 1;
for n_pass = 1:1:20 %Iterate through number of passengers
    for ban_length = 10/12*ftTOm : 0.1 : 5*ftTOm %WHERE DOES BANNER LENGTH GET FACTORED IN!!!!!
        for v_c = 1:0.05:100 %Iterate through cruise velocities
            for b = 1*ftTOm:0.1:5*ftTOm %Iterate through wingspan
                n_lugg = n_pass;
                m_comp = n_pass*m_pass + n_lugg*m_lugg;
                m_tot = m_motor + m_bat + m_struct + m_comp + m_wing;
                p = aerodynamics(ae, m_tot, v_c, Cf_ban, ban_length, c);
                F_D = 0.5*c.rho*v_c^2 * ae.S*p.CD;
                P_aero = v_c*F_D;
                n_laps = 0;
        %         disp(P_aero);
        %         disp(P_mech);
                takeoff()
                score = 0;
                if per_diff(P_aero, P_mech) <= 0.0005
                    configs{end+1} = {P_aero, P_mech, n_pass, v_c};
                end
            end
        end
    end
end
%Find where mechanical and aerodynamic power match, ie. a successful match
% diff = abs(P_aero_ - P_mech_);
% [val, ind] = min(diff(:));
% [row, col] = ind2sub([size(P_aero_,1), size(P_aero_,2)], ind);
% Power = P_aero_(row, col);

function ae = aircraft(b, c_root, lam, e, CD0, al)
    %Calculates properties of an aircraft
    ae.S = 0.5*b*c_root*(1+lam); %Planform area of tapered wing
    ae.AR = b^2/ae.S;
    ae.b = b;
    ae.c_root = c_root;
    ae.lam = lam;
    ae.e = e;
    ae.CD0 = CD0;
    ae.al = al;
end

function p = aerodynamics (ae, m, v, Cf, ban_length, c)
    %b - full wing span
    %c_root - root chord
    %lam - taper ratio = tip/root chord
    %e - Oswald efficiency
    %rho - air density
    %m - total aircraft mass
    %v - velocity
    %CD0 - zero-lift coefficient of drag (determined experimentally)
    p.K = 1/(pi*ae.e*ae.AR);
    p.Q = 0.5*c.rho*v^2; %Dynamic pressure
    p.CL = m*c.g/(p.Q*ae.S); %CL = L/(Q*S)
    p.CD = ae.CD0 + p.K*p.CL^2 + Cf;
end

function score = calcScore (n_laps, b_length)
    M3 = 2 + n_laps*b_length;
    score = M3;
end

function diff = per_diff(val1, val2)
    diff = abs(val1 - val2)/val1;
end

function t_to = takeoff(ae, m)
    z0 = [0, 0]; %initial velocity and speed
    t = linspace(1, 100, 1000);
    myFun = @(z,t, ae, m) odes(z, t, ae, m);
    [t_array, z_array] = ode45(myFun, z0, t, odeset('RelTol', 1e-10, 'AbsTol', 1e-10));
    
    function zdot = odes (z, t, ae, m)
        x = z(1); v = z(2);
        p = aerodynamics(ae, m, v, Cf, ban_length, c);
        xdot = v;
        vdot = 1/m * (T*cos(ae.al) - 0.5*rho*v^2*ae.S*p.CD - c.mu*(m*g - 0.5*c.rho*v^2*ae.S*p.CD));
        zdot = [xdot, vdot];
    end
end