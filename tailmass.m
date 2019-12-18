%TAIL WEIGHT SCRIPT (BASIC) 
%This script is a simple calculation script for tail mass. (may be modified
%for future use) 

t = input('thickness of tail = ');
b = input('span of tail = ');
c = input('chord of tail = '); 
d = input('density of tail material = ');

% MAKE SURE UNITS MATCH!

v = t*b*c; %Volume of tail

M = d*v; %Mass of tail

fprintf('Mass is %f\n', M)