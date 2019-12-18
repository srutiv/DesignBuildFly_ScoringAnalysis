function aircraft=generate_standard_aircraft(plane_name)
% generate standard aircraft takes in a string, which is to be the
% aircraft's name, and generates an Aircraft object with properties that
% are a normal distribution about a pre-determined mean aircraft

%% determine wing values 
aileron=Surface(.75,-1,.1,.9,'aileron');
b=normrnd(1,.125)*10;
S=normrnd(1,.125)*8;
taper=rand;
sweep=rand*pi/8;
airfoil='mh114.dat';
dihedral=rand*pi/20;
surfaces=[aileron];
coord=[0 0 0];
name='wing';

wing=Main_Wing(b,S,taper,sweep,airfoil,airfoil,dihedral,surfaces,coord,name)

%% Horizontal Tail
elevator=Surface(.75,1,.1,.9,'elevator');
bwing=b;
Swing=S;
b=normrnd(1,.125)*4;
S=normrnd(1,.125)*2;
taper=rand;
sweep=rand*pi/4;
airfoil=[];
dihedral=rand*pi/20;
surfaces=elevator;
croot=wing.chord(0)
coord=[2*croot+croot*normrnd(0,.25) 0 0];
name='htail';

htail=Hor_Stab(b,S,taper,sweep,airfoil,dihedral,surfaces,coord,name);
%% Vertical tail
rudder=Surface(.75,1,.1,.9,'rudder')
b=normrnd(1,.125)*1.5;
S=normrnd(1,.125)*1.5;
taper=rand;
sweep=rand*pi/4;
airfoil=[];
surfaces=[rudder];
name='vtail';
vtail=Ver_Stab(b,S,taper,sweep,airfoil,surfaces,coord,name);
%% Generate Plane
M=Swing*.5+htail.s*.5+vtail.s*.5
CG=[Swing/bwing*.5 0 0];
IXX=bwing^2*M;
IYY=htail.coord(1)^2*M+croot^2*M;
IZZ=IXX+IYY;
%% initialize test parameters 
param.d=1.225;
param.v=10;
param.g=9.6;
param.use_alpha=0;
param.alpha=0;
param.beta=0;
aircraft=Aircraft([wing,wing],htail,vtail,[],CG,M,[IXX IYY IZZ],plane_name,param);

end