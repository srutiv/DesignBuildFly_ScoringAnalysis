% Michael Zakoworotny
% DBF
% Preliminary Scoring Testing
clear all; clc;
%Conversions
ftTOm = 0.3048; % ft to m conversion
ozTOkg = 0.0283495; %oz to kg covnersion

%First approximations of geometric values
% b = 5 * ftTOm; %m, span
t_3 = 10 * 60; %seconds, M3 duration
d_takeoff = 20 * ftTOm; %m, takeoff distance
d_lap = 3570.79 * ftTOm; %m, lap distance
c_root = 1 * ftTOm; %root chord
lam = 0.9; %taper ratio
e = 0.8; %Oswald efficiency, for "other aircraft" according to Raymer
CD0 = 0.02; %Zero-lift coefficient of drag, for a "clean propeller aircraft" according to Raymer
CL_max = 1.20; %Maximum coefficient of lift - ie. CL at stall angle
al = 0; %Angle of attack
m_pass = 4 * ozTOkg; %mass of passenger
m_lugg = 1 * ozTOkg; %mass of luggage
m_struct = 1;
m_wing = 0.5;
Cf_ban = 7.72e-3; %Cf = Cd for objects with no pressure drag

%Lithium Ion Data, https://www.claytonpower.com/products/lithium-ion-batteries/specifications/
q = 100 * 3600; %Amp*hours -> Amp*second
V = 12.6; % Volts
I = 63; %Continuous discharge current
eta = 0.7/2; %efficiency of motor
P_motor = I*V; %Power of motor
% batt_dens = 10 * P_motor/1000*3600; %kg/kWh
m_bat = 178.03/1000; %kg mass of battery, very rought estimate
m_motor = 3.10 * ozTOkg; %cobra motor mass, https://www.cobramotorsusa.com/motors-2221-16.html
P_mech_max = eta*P_motor;

%Constants
c.rho = 1.14; %kg/m^3, calculated density in Wichita using altitude
c.visc = 1.765e-5;
c.mu = 0.02;
c.g = 9.81;

% n_pass = 1:1:20;
% v_c = 1:0.5:100;
configs = {};
i = 1;
% for n_pass = 1:1:20 %Iterate through number of passengers
% for P_mech = 0 : 0.1 : P_mech_max
    for ban_length = 10/12*ftTOm : 1*ftTOm : 5*ftTOm %WHERE DOES BANNER LENGTH GET FACTORED IN!!!!!
%         for P = 0 : 0.01 : I*V
        for v_c = 1:5:100 %Iterate through cruise velocities
            for b = 1*ftTOm:1*ftTOm:5*ftTOm %Iterate through wingspan
%                 n_lugg = n_pass;
%                 m_comp = n_pass*m_pass + n_lugg*m_lugg;
                m_tot = m_motor + m_bat + m_struct + m_wing;% + m_comp;
                ap = aircraft(b, c_root, lam, e, CD0, al, CL_max, m_tot, c); %Store these parameters in a struct
                pr = prop(ap,V,I,eta);
                [t_to, x_to, v_to] = takeoff(ap, pr, m_tot, Cf_ban, ban_length, c, d_takeoff);
                aer = aerodynamics(ap, m_tot, v_c, Cf_ban, ban_length, c, 0);
                
%                 F_D = 0.5*c.rho*v_c^2 * ap.S*aer.CD;
%                 P_aero = v_c*F_D;
        %         disp(P_aero);
        %         disp(P_mech);
                aer.CL;
                
                if x_to <= d_takeoff
                    if 0.5*c.rho*v_c^2*ap.S*aer.CL >= m_tot*c.g %Check if lift makes sense
                        if 0.5*c.rho*v_c^2*ap.S*aer.CD <= pr.T
                            L = 0.5*c.rho*v_c^2*ap.S*aer.CL;
                            W = m_tot*c.g;
                            t_flight = t_3 - t_to;
                            nlaps = v_c*t_flight / d_lap;
                            sc = score(t_3, t_to, nlaps, ban_length);
            %                 if per_diff(P_aero, P_mech) <= 0.0005
        %                         configs{end+1} = {P_aero, P_mech, v_c, sc};
            %                 end
                            fprintf("Score: %.f, Vc: %.2f, nlaps: %.2f, banL: %.2f, v_to: %.2f\n", sc, v_c, nlaps, ban_length, v_to);
                        else
                            disp("nope");
                        end
                    end
                end
            end
        end
    end
% end
% end
%Find where mechanical and aerodynamic power match, ie. a successful match
% diff = abs(P_aero_ - P_mech_);
% [val, ind] = min(diff(:));
% [row, col] = ind2sub([size(P_aero_,1), size(P_aero_,2)], ind);
% Power = P_aero_(row, col);

function ap = aircraft(b, c_root, lam, e, CD0, al, CL_max, m, c)
    %Calculates properties of an aircraft
    %b - full wing span
    %c_root - root chord
    %lam - taper ratio = tip/root chord
    %e - Oswald efficiency
    %CD0 - zero-lift coefficient of drag (determined experimentally)
    %al - angle of attack
    %CL_max - max CL, used for takeoff
    %m_tot - total mass of aircraft
    %c - struct containing constants
    ap.S = 0.5*b*c_root*(1+lam); %Planform area of tapered wing
    ap.AR = b^2/ap.S;
    ap.b = b;
    ap.c_root = c_root;
    ap.lam = lam;
    ap.e = e;
    ap.CD0 = CD0;
    ap.al = al;
    ap.CL_max = CL_max;
    ap.vs = sqrt(m*c.g / (0.5*c.rho*ap.S*ap.CL_max));
end

function aer = aerodynamics (ap, m, v, Cf, ban_length, c, ground)
    %ap - struct containing aircraft properties
    %m - total mass of plane
    %v - current velocity
    %Cf - coefficient of friction of banner
    %ban_length - length of banner
    %c - struct containing constants
    %ground - boolean specifying current phase: 1 for ground, 0 for cruise
    aer.Vs = sqrt(m*c.g / (0.5*c.rho*ap.S*ap.CL_max));
%     aer.vc = (2*P/(c.rho*ap.S))^(1/3)
    aer.K = 1/(pi*ap.e*ap.AR);
    aer.Q = 0.5*c.rho*v^2; %Dynamic pressure
    if ground %Use CLmax before takeoff
        aer.CL = ap.CL_max;
        aer.CD = ap.CD0 + aer.K*ap.CL_max^2 + Cf*ban_length;
    else %Use calculated CL during cruise
        aer.CL = m*c.g/(aer.Q*ap.S); %CL = L/(Q*S)
        aer.CD = ap.CD0 + aer.K*aer.CL^2 + Cf*ban_length;
    end
end

function pr = prop(ap,V,I,eta)
    P_max = V*I;
    P_mech = eta*P_max;
    pr.T = P_mech/ap.vs;
end

function diff = per_diff(val1, val2)
    diff = abs(val1 - val2)/val1;
end

function [t_takeoff, x_takeoff, v_takeoff] = takeoff(ap, pr, m, Cf, ban_length, c, d_takeoff)
    z0 = [0, 0]; %initial position and speed
    t = [0, 100];
    myFun = @(t,z) odes(t, z, ap, pr, m, Cf, ban_length, c);
    [t_, z_] = ode45(myFun, t, z0, odeset('RelTol', 1e-10, 'AbsTol', 1e-10));
    x_ = z_(:,1); v_ = z_(:,2);
    v_takeoff = 1.3*ap.vs;
    t_takeoff = interp1(v_, t_, v_takeoff); %Interpolates to find the time at which it reaches takeoff
    x_takeoff = interp1(t_, x_, t_takeoff);
    
    function zdot = odes (t, z, ap, pr, m, Cf, ban_length, c)
        x = z(1); v = z(2);
        aer = aerodynamics(ap, m, v, Cf, ban_length, c, 1);
        xdot = v;
        vdot = 1/m * (pr.T*cos(ap.al) - 0.5*c.rho*v^2*ap.S*aer.CD - c.mu*(m*c.g - 0.5*c.rho*v^2*ap.S*aer.CL));
        zdot = [xdot; vdot];
    end
end

function sc = score(t_3, t_to, nlaps, ban_length)
    
    sc = nlaps * ban_length;
end