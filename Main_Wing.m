classdef Main_Wing
    properties
        b %full_span
        sb % full straight span
        s %area
        taper % taper tatio, c_tip/c_root
        sweep % sweep angle in radians (0 is no sweep)
        airfoil % airfoil file name (i.e. 'airfiol.dat')
        dihedral % dihedral angle in radians (0 is no dihedral)
        angle % incidence angle in degrees
        surfaces % array of surface class objects
        coord % coordinates of the center of the leading edge i.e. [X,Y,Z]
        name
        %%%% Y offset should be 0, unless two main wings are created to
        %%%% make wingtips, or a wing that looks like \___/ (i.e pelikan
        %%%% tail)
    end
    
    methods
%% CONSTRUCTOR
        function obj=Main_Wing(B,SB,S,Taper,Sweep,Airfoil,Dihedral,Incidence,Surfaces,Coord,Name) %% class constructor
            obj.b=B;
            obj.sb=SB;
            obj.s=S;
            obj.taper=Taper;
            obj.sweep=Sweep;
            obj.airfoil=Airfoil;
            obj.dihedral=Dihedral;
            obj.angle=Incidence;
            obj.surfaces=Surfaces;
            obj.coord=Coord;
            obj.name=Name;
        end
        
%% CHORD at Y (AVL axis -- y goes along span)
% Chord as a function of span
% y = 0 is the center of the plane like in AVL
function cy=chord(obj,y)
    %cy=2*obj.s/(1+obj.taper)/obj.b*(1-(2*(1-obj.taper))/obj.b*y); %for
    %normal planes
    
    % for twin boom plane
    croot = 2*obj.s/(1+obj.taper)/obj.b*(1-(2*(1-obj.taper))/obj.b*0);
    if (y < obj.sb/2)
        cy = croot; %root chord
    else
        s = obj.s-(croot*obj.sb); y = y - obj.sb/2; b = obj.b - obj.sb;
        cy=2*s/(1+obj.taper)/b*(1-(2*(1-obj.taper))/b*y);
    end
end
      
%% MEAN CHORD

function MAC=mean_chord(obj)
    l=obj.taper;
    b=obj.b;
    s=obj.s;
    y=b/2;;
    MAC=(2/s)*(2*s/b/(1+l))^2*(y-2*(1-l)*y^2/b+4/3*(1-l)^2*y^3/b^2);
end

%% PLOT WIREFRAME WING
        function plot_wing(obj)
            %% plots the wing as 6 points on a 3 axis graph. Sets axis to equal
            p1x=obj.coord(1);
            p1y=obj.coord(2);
            p1z=obj.coord(3);
            
            croot=chord(obj,0);
            ctip=chord(obj,obj.b/2);
            
            p2x=p1x+croot;
            p2y=p1y;
            p2z=p1z;
            
            p3x=(p1x+croot*.25)+tan(obj.sweep)*obj.b/2-ctip*.25;
            p3y=p1y+obj.b/2;
            p3z=p1z+obj.b/2*sin(obj.dihedral);
            
            p4x=p3x+ctip;
            p4y=p3y;
            p4z=p3z;
            
            X=[p1x,p2x,p4x,p3x,p1x,p3x,p4x,p2x];
            Y=[p1y,p2y,p4y,p3y,p1y,-p3y,-p4y,p2y];
            Z=[p1z,p2z,p4z,p3z,p1z,p3z,p4z,p2z];
            hold on
            plot3(X,Y,Z,'k')
            
            if ~isempty(obj.surfaces)
                for n=1:length(obj.surfaces)
                    surface=obj.surfaces(n);
                    chinge=surface.chinge;
                    bi=surface.bi;
                    be=surface.be;
                    ps1x=p1x+croot*.25+tan(obj.sweep)*obj.b/2*bi+.75*chord(obj,obj.b*bi/2);
                    ps1y=p1y+obj.b/2*bi;
                    ps1z=p1z+obj.b/2*sin(obj.dihedral)*bi;
                    
                    ps2x=p1x+croot*.25+tan(obj.sweep)*obj.b/2*be+.75*chord(obj,obj.b*be/2);
                    ps2y=p1y+obj.b/2*be;
                    ps2z=p1z+obj.b/2*sin(obj.dihedral)*be;
                    
                    ps3x=p1x+croot*.25+tan(obj.sweep)*obj.b/2*be+(chinge-.25)*chord(obj,obj.b*be/2);
                    ps3y=p1y+obj.b/2*be;
                    ps3z=p1z+obj.b/2*sin(obj.dihedral)*be;
                    
                    ps4x=p1x+croot*.25+tan(obj.sweep)*obj.b/2*bi+(chinge-.25)*chord(obj,obj.b*bi/2);
                    ps4y=p1y+obj.b/2*bi;
                    ps4z=p1z+obj.b/2*sin(obj.dihedral)*bi;
                    X=[ps1x,ps2x,ps3x,ps4x,ps1x];
                    Y=[ps1y,ps2y,ps3y,ps4y,ps1y];
                    Z=[ps1z,ps2z,ps3z,ps4z,ps1z];
                    plot3(X,Y,Z,'r')
                    plot3(X,-Y,Z,'r')
                end
            end

        end
%% BUILD AVL FILE
        function build_surface(obj,fid)
            surfaces=obj.surfaces;
            fprintf(fid,'#========================================== \n');
            fprintf(fid,'SURFACE \n');
            fprintf(fid,'%s \n',obj.name);
            fprintf(fid,'#Nchordwise Cspace Nspanwise Sspace \n');
            fprintf(fid,'12 1.0 16 1.0 \n');
            fprintf(fid,'# \n');
            fprintf(fid,'YDUPLICATE \n');
            fprintf(fid,'0.0 \n');
            fprintf(fid,'TRANSLATE\n');
            fprintf(fid,'%f %f %f\n',obj.coord(1),obj.coord(2),obj.coord(3));
            fprintf(fid,'ANGLE \n');
            fprintf(fid,'%f \n',obj.angle);
            % section calculations
            n_sections=2+2*length(surfaces);
            section_array=[0,1];
            for n=1:length(surfaces);
                section_array(n+1,1)=surfaces(n).bi;
                section_array(n+1,2)=surfaces(n).be;
            end
            
            croot=chord(obj,0);
            
            sections=unique(section_array);
            for n=1:length(sections)
                fprintf(fid,'#------------------------------------------ \n');
                fprintf(fid,'SECTION \n');
                fprintf(fid,'#Xle Yle Zle Chord Ainc Nspanwise Sspace \n');
                b_frac=sections(n);
                c=chord(obj,obj.b*b_frac/2);
                xn=croot*.25+tan(obj.sweep)*obj.b*b_frac/2-c*.25;
                yn=b_frac*obj.b/2;
                zn=b_frac*obj.b/2*sin(obj.dihedral);
                fprintf(fid,'%f %f %f %f 0.0 0 0 \n',xn,yn,zn,c);
                fprintf(fid,'AFILE \n');
                fprintf(fid,'%s \n',obj.airfoil);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for o=1:length(surfaces) 
                    disp(o)
                    if surfaces(o).bi==b_frac || surfaces(o).be==b_frac
                        fprintf(fid,'#Cname Cgain Xhinge HingeVec   SgnDup\n');
                        fprintf(fid,'CONTROL \n');
                        fprintf(fid,'%s 1 %f 0.0 0.0 0.0 %f\n',surfaces(o).name,surfaces(o).chinge,surfaces(o).sign);
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
        end
%% PLOT AIRFOIL SHAPE     
        function plot_airfoil(obj)
            % plots the wing as 6 points on a 3 axis graph. Sets axis to equal
            a=importdata(obj.airfoil);
            a_array=a.data;
            p1x=obj.coord(1);
            p1y=obj.coord(2);
            p1z=obj.coord(3);
            croot=chord(obj,0);
            hold on
            N=obj.b*10;
            for n=0:N
                b_frac=n/N;
                c=chord(obj,obj.b*b_frac/2);
                xn=p1x+croot*.25+tan(obj.sweep)*obj.b*b_frac/2-c*.25;
                yn=p1y+b_frac*obj.b/2;
                zn=p1z+b_frac*obj.b/2*sin(obj.dihedral);
                plot3(xn+c*a_array(:,1),yn+a_array.*0,zn+c*a_array(:,2),'k')
                plot3(xn+c*a_array(:,1),-yn-a_array.*0,zn+c*a_array(:,2),'k')
                if mod(n,2)==0
                    if ~isempty(obj.surfaces)
                        for o=1:length(obj.surfaces)
                            surface=obj.surfaces(o);
                            chinge=surface.chinge;
                            bi=surface.bi;
                            be=surface.be;
                            new_array=[0,0];
                            if (n/N >= bi) && (n/N <= be)
                                p=n/N;
                                q=1;
                                for i=1:length(a_array(:,2))
                                    if a_array(i,1) >= chinge
                                        new_array(q,:)=a_array(i,:);
                                        q=q+1;
                                    end
                                end
                                plot3(xn+c*new_array(:,1),yn+new_array.*0,zn+c*new_array(:,2),'r')
                                plot3(xn+c*new_array(:,1),-yn+new_array.*0,zn+c*new_array(:,2),'r')
                                clear new_array
                            end
                        end
                    end
                end
                
                axis equal
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAKE TRAILING EDGE HORIZONTAL
% function obj=hor_trail(obj)
%     obj.sweep=asin(2/obj.b*(obj.croot*(1-obj.taper)-.25*obj.croot+obj.croot*obj.taper*.25));
% end
    end
end
