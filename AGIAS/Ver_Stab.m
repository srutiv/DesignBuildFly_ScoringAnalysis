classdef Ver_Stab
    %% PROPERTIES
    properties
        b %span
        s %root chord
        taper % taper tatio, c_tip/c_root
        sweep % sweep angle in radians (0 is no sweep)
        airfoil % airfoil file name (i.e. 'airfiol.dat')
        surfaces % array of surface class objects
        coord % coordinates of the center of the leading edge i.e. [X,Y,Z]
        name
    end
    %% METHODS
    methods
        function obj=Ver_Stab(B,S,Taper,Sweep,Airfoil,Surfaces,Coord,Name) %% class constructor
            obj.b=B;
            obj.s=S;
            obj.taper=Taper;
            obj.sweep=Sweep;
            obj.airfoil=Airfoil;
            obj.surfaces=Surfaces;
            obj.coord=Coord;
            obj.name=Name;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% CHORD at Y
        function cz=chord(obj,z)
            cz=obj.s*2/(1+obj.taper)/obj.b*(1-((1-obj.taper))/obj.b*z);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_wing(obj)
            %% plots the wing as 6 points on a 3 axis graph. Sets axis to equal
            p1x=obj.coord(1);
            p1y=obj.coord(2);
            p1z=obj.coord(3);
            
            croot=chord(obj,0);
            ctip=chord(obj,obj.b);
            
            p2x=p1x+croot;
            p2y=p1y;
            p2z=p1z;
            
            p3x=p1x+croot*.25+tan(obj.sweep)*obj.b-ctip*.25;
            p3y=p1y;
            p3z=p1z+obj.b;
            
            p4x=p3x+ctip;
            p4y=p1y;
            p4z=p3z;
            
            X=[p1x,p2x,p4x,p3x,p1x];
            Y=[p1y,p2y,p4y,p3y,p1y];
            Z=[p1z,p2z,p4z,p3z,p1z];
            hold on
            plot3(X,Y,Z,'k')
            X=[p1x,p2x,p4x,p3x];
            Y=[p1y,p2y,p4y,p3y];
            Z=[p1z,p2z,p4z,p3z];
            color=.8;
            fill3(X,Y,Z,[color color color])
            fill3(X,-Y,Z,[color color color])
            
            if ~isempty(obj.surfaces)
                for n=1:length(obj.surfaces)
                    surface=obj.surfaces(n);
                    chinge=surface.chinge;
                    bi=surface.bi;
                    be=surface.be;
                    ps1x=p1x+croot*.25+tan(obj.sweep)*obj.b*bi+.75*chord(obj,obj.b*bi);
                    ps1y=p1y;
                    ps1z=p1z+obj.b*bi;
                    
                    ps2x=p1x+croot*.25+tan(obj.sweep)*obj.b*be+.75*chord(obj,obj.b*be);
                    ps2y=p1y;
                    ps2z=p1z+obj.b*be;
                    
                    ps3x=p1x+croot*.25+tan(obj.sweep)*obj.b*be+(chinge-.25)*chord(obj,obj.b*be);
                    ps3y=p1y;
                    ps3z=p1z+obj.b*be;
                    
                    ps4x=p1x+croot*.25+tan(obj.sweep)*obj.b*bi+(chinge-.25)*chord(obj,obj.b*bi);
                    ps4y=p1y;
                    ps4z=p1z+obj.b*bi;
                    X=[ps1x,ps2x,ps3x,ps4x,ps1x];
                    Y=[ps1y,ps2y,ps3y,ps4y,ps1y];
                    Z=[ps1z,ps2z,ps3z,ps4z,ps1z];
                    plot3(X,Y,Z,'r')
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_airfoil(obj)
            %% plots the wing as 6 points on a 3 axis graph. Sets axis to equal
            a=importdata(obj.airfoil);
            a_array=a.data;
            p1x=obj.coord(1);
            p1y=obj.coord(2);
            p1z=obj.coord(3);
            hold on
            N=obj.b*20;
            new_array=[0,0];
            tic
            for n=0:N
                b_frac=n/N;
                c=chord(obj,obj.b*b_frac);
                xn=p1x+croot*.25+tan(obj.sweep)*obj.b*b_frac/2-c*.25;
                yn=p1y+b_frac*obj.b/2;
                zn=p1z+n/N*obj.b;
                plot3(xn+c*a_array(:,1),yn+a_array(:,2),zn+a_array(:,2).*0,'k')
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
                                plot3(xn+c*new_array(:,1),yn+new_array(:,2),zn+new_array(:,2).*0,'r')
                                
                                clear new_array
                            end
                        end
                    end
                end
                
                axis equal
            end
            toc
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function build_surface(obj,fid)
            surfaces=obj.surfaces;
            fprintf(fid,'#========================================== \n');
            fprintf(fid,'SURFACE \n');
            fprintf(fid,'%s \n',obj.name);
            fprintf(fid,'#Nchordwise Cspace Nspanwise Sspace \n');
            fprintf(fid,'12 1.0 16 1.0 \n');
            fprintf(fid,'TRANSLATE\n');
            fprintf(fid,'%f %f %f\n',obj.coord(1),obj.coord(2),obj.coord(3));
            fprintf(fid,'YDUPLICATE \n');
            fprintf(fid,'0.0 \n');
            n_sections=2+2*length(surfaces);
            section_array=[0,1];
            for n=1:length(surfaces)
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
                c=chord(obj,obj.b*b_frac);
                xn=croot*.25+tan(obj.sweep)*obj.b*b_frac-c*.25;
                yn=0;
                zn=b_frac*obj.b;
                fprintf(fid,'%f %f %f %f 0.0 0 0 \n',xn,yn,zn,c);
                if ~isempty(obj.airfoil)
                    fprintf(fid,'AFILE \n');
                    fprintf(fid,'%s \n',obj.airfoil);
                end
                for o=1:length(surfaces)
                    disp(o)
                    if surfaces(o).bi==b_frac || surfaces(o).be==b_frac
                        fprintf(fid,'#Cname Cgain Xhinge HingeVec   SgnDup\n');
                        fprintf(fid,'CONTROL \n');
                        fprintf(fid,'%s 1 %f 0.0 0.0 0.0 %f\n',surfaces(o).name,surfaces(o).chinge,surfaces(o).sign);
                    end
                end
                
            end
        end
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function obj=flush_trail(obj)
        %     obj.sweep=atan(1/obj.b*(obj.croot*(1-obj.taper)-.25*obj.croot+obj.croot*obj.taper*.25));
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end