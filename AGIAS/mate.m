function pout=mate(p1,p2,add_num)
% mate takes in two Aircraft objects and an integer, and outputs an
% aircraft whose parameters are a randomly generated weighted sum of the
% two input aircraft.


Wing1=p1.wings;
Htail1=p1.htails;
Vtail1=p1.vtails;
Wing2=p2.wings;
Htail2=p2.htails;
Vtail2=p2.vtails;
Wing=Wing1;
Htail=Htail1;
Vtail=Vtail1;


Wing.b=(Wing1.b+Wing2.b)/2*normrnd(1,.25);
Wing.s=(Wing1.s+Wing2.s)/2*normrnd(1,.25);
Wing.taper=(Wing1.taper+Wing2.taper)/2*normrnd(1,.25);
Wing.sweep=(Wing1.sweep+Wing2.sweep)/2*normrnd(1,.25);
Wing.dihedral=(Wing1.dihedral+Wing2.dihedral)/2*normrnd(1,.25);

Htail.b=(Htail1.b+Htail2.b)/2*normrnd(1,.25);
Htail.s=(Htail1.s+Htail2.s)/2*normrnd(1,.25);
Htail.taper=(Htail1.taper+Htail2.taper)/2*normrnd(1,.25);
Htail.sweep=(Htail1.sweep+Htail2.sweep)/2*normrnd(1,.25);
Htail.coord(1)=(Htail1.coord(1)+Htail2.coord(1))/2*normrnd(1,.25);

Vtail.b=(Vtail1.b+Vtail2.b)/2*normrnd(1,.25);
Vtail.s=(Vtail1.s+Vtail2.s)/2*normrnd(1,.25);
Vtail.taper=(Vtail1.taper+Vtail2.taper)/2*normrnd(1,.25);
Vtail.sweep=(Vtail1.sweep+Vtail2.sweep)/2*normrnd(1,.25);
Vtail.coord(1)=Htail.coord(1);
bwing=Wing.b;
cwing=Wing.mean_chord;

M=Wing.s*.5+Htail.s*.125+Vtail.s*.125;
CG=[cwing*.5 0 0];
IXX=bwing^2*M;
IYY=Htail.coord(1)^2*M;
IZZ=IXX+IYY;
Name=sprintf('%s_%i',p1.name,add_num);
Inert=[IXX IYY IZZ];

pout=Aircraft(Wing,Htail,Vtail,[],CG,M,Inert,Name);

end