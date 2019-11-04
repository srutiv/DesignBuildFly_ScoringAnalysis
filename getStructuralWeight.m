function mass_struct = getStructuralWeight(d, l, b, c)
% Inputs:
% d: Nominal diameter of fuselage (for convservative estimate, assume
% fuselage has a square cross-section
% l: length of fuselage
% b: wing span
% c: mean aerodynamic chord

% REFERENCE POINTS
% Iteration 1 Plane
% d = mean([108, 120]) = 114mm
% l = 1050mm
% b = 1500 mm
% c = 480 mm
% Fuselage mass = 141.02 g
% Wing mass = 316.66 g

% Example:
% getStructuralWeight(114, 1050, 1500, 480) returns the weight in grams of
% the main structural components of our iteration 1 aircraft

% TODO:
% - Add more data points (use Anthony's sheet with previous competition
% planes)
% - Refine stastical distribution

scale = 3000; %scale mass of plane (implemented with probability distribution

refFus = 141.02;
refWing = 316.66;
refD = mean([108,120]);
refL = 1050;
refB = 1500;
refC = 480;

fusMass = binornd(scale, refFus/scale) * (d/refD)^2 * (l/refL);
wingMass = binornd(scale, refWing/scale) * (b/refB) * (c/refC);

mass_struct = (fusMass + wingMass)/1000;
