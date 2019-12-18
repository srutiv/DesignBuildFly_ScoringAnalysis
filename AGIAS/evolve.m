function planeout=evolve(planein,add_num)
% takes inputs of an Aircraft object, and an integer to be added to the
% aircrafts name. Outputs a new aircraft with characteristics that are
% randomly varied with a normal distribution.
Wing=planein.wings;
Htail=planein.htails;
Vtail=planein.vtails;

Wing.b=Wing.b*normrnd(1,.125);
Wing.s=Wing.s*normrnd(1,.125);
Wing.taper=Wing.taper+normrnd(0,.125);
if Wing.taper < 0 
   Wing.taper = 0 ;
end
if Wing.taper > 1
    Wing.taper=1;
end
Wing.sweep=Wing.sweep*normrnd(1,.125);
Wing.dihedral=Wing.dihedral*normrnd(1,.125);

Htail.b=Htail.b*normrnd(1,.125);
Htail.s=Htail.s*normrnd(1,.125);
Htail.taper=Htail.taper+normrnd(0,.125);
Htail.dihedral=Htail.dihedral*normrnd(1,.125);
if Htail.taper < 0 
   Htail.taper = 0 ;
end
Htail.sweep=Htail.sweep*.75+Htail.sweep*rand*.5;
x=normrnd(0,.125)*Htail.coord(1)+Htail.coord(1);
Htail.coord(1)=x;

Vtail.b=Vtail.b*normrnd(1,.125);
Vtail.s=Vtail.s*normrnd(1,.125);
Vtail.taper=Vtail.taper+normrnd(0,.125);
if Vtail.taper < 0 
   Vtail.taper = 0 ;
end
Vtail.sweep=Vtail.sweep*normrnd(1,.125);
Vtail.coord(1)=x;

bwing=Wing.b;
cwing=Wing.mean_chord;

M=Wing.s*.5+Htail.s*.125+Vtail.s*.125;
CG=[cwing*.5 0 0];
IXX=bwing^2*M;
IYY=Htail.coord(1)^2*M+cwing^2*M;
IZZ=IXX+IYY;
Name=sprintf('%s_%i',planein.name,add_num);
Inert=[IXX IYY IZZ];

planeout=Aircraft(Wing,Htail,Vtail,[],CG,M,Inert,Name);
end