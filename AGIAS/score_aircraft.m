function score=score_aircraft(score_matrix,plane)

%% Extract Aircraft Parameters 
st=plane.st;
run=plane.run;
CL=run.CLtot;
CD=run.CDtot;
CL_CD=CL/CD;
e=run.e;
Clb=st.Clb; %% IMPORTANT roll moment relative to side slip
Clp=st.Clp; %% IMPORTANT roll moment relative to roll rate
Clr=st.Clr;
Cnb=st.Cnb; %% IMPORTANT yaw moment side slip
Cnp=st.Cnp; 
Cnr=st.Cnr; %% IMPORTANT yaw moment yaw rate
Cma=st.Cma; %% IMPORTANT pitch moment angle of attack
Cmq=st.Cmq;
wings=plane.wings(1);
AR=wings.b^2/wings.s;
L=CL*wings.s*.5*1.225*25^2;
roll=plane.eig.roll(1);
dutch=plane.eig.dutch(1);
short=plane.eig.short(1);
spiral=plane.eig.spiral(1);
phugoid=plane.eig.phugoid(1);

aircraft_values=[CL;CD;CL_CD;e;Clb;Clp;Clr;Cnb;Cnp;Cnr;Cma;Cmq;AR;L;roll;dutch;short;spiral;phugoid];
aircraft_scores=abs(aircraft_values-score_matrix(:,1)).*score_matrix(:,2);
score=sum(aircraft_scores);

end