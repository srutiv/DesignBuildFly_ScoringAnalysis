classdef Aircraft
    properties
        wings %% array of Main_Wing objects
        htails %% array of Hor_Tail objects
        vtails %% array of Ver_Tail objects
        bodies %% array of Body Objects
        cg %% Mass, CG, and moments of inertia [X Y Z]
        m %% mass of plane
        inert %% IXX IYY IZZ
        name %% name of the aircraft
        st  %% struct that gives the stability derivatives
        run %% struct that gives run information 
        eig %% roots of plane's stability matrix. From top to bottom: Roll, Dutch Roll, Short, Slip, Phugoid
        param
        %%col 1 is real part, col 2 is positive imaginary root.
    end
    methods
        %% CONSTRUCTOR
        function obj=Aircraft(Wings,Htails,Vtails,Bodies,CG,M,Inert,Name,param)
            obj.wings=Wings;
            obj.htails=Htails;
            obj.vtails=Vtails;
            obj.bodies=Bodies;
            obj.cg=CG;
            obj.m=M;
            obj.inert=Inert;
            obj.name=Name;
            obj.st=[];
            obj.run=[];
            obj.eig.roll=[0,0];
            obj.eig.dutch=[0,0];
            obj.eig.short=[0,0];
            obj.eig.spiral=[0,0];
            obj.eig.phugoid=[0,0];
            obj.param=param;
        end
        %% Plot Airfoils
        function plot_all_airfoil(obj)
            if ~isempty(obj.wings)
                for n=1:length(obj.wings)
                    wing=obj.wings(n);
                    wing.plot_airfoil;
                end
            end
            if ~isempty(obj.htails)
                for n=1:length(obj.wings)
                    tail=obj.htails(n);
                    tail.plot_airfoil;
                end
            end
            if ~isempty(obj.vtails)
                for n=1:length(obj.wings)
                    tail=obj.vtails(n);
                    tail.plot_airfoil;
                end
            end
            if ~isempty(obj.bodies)
                for n=1:length(obj.bodies)
                    body=obj.bodies(n);
                    body.plot_body
                end
            end
        end
        %% Plot Wireframe
        function plot_all(obj)
            if ~isempty(obj.wings)
                for n=1:length(obj.wings)
                    wing=obj.wings(n);
                    wing.plot_wing;
                end
            end
            if ~isempty(obj.htails)
                for n=1:length(obj.htails)
                    tail=obj.htails(n);
                    tail.plot_wing;
                end
            end
            if ~isempty(obj.vtails)
                for n=1:length(obj.vtails)
                    tail=obj.vtails(n);
                    tail.plot_wing;
                end
            end
            if ~isempty(obj.bodies)
                for n=1:length(obj.bodies)
                    body=obj.bodies(n);
                    body.plot_body
                end
            end
            if ~isempty(obj.cg)
                plot3(obj.cg(1),obj.cg(2),obj.cg(3),'*r')
            end
            if ~isempty(obj.st)
                plot3(obj.st.NP,0,0,'*b')
            end
            axis equal
        end
        %% BUILD AVL FILE
        function build_file(obj)
            string=['.\AVL_FILES\',obj.name,'.avl'];
            fid=fopen(string,'wt');
            fid=fclose(fid);
            wing=obj.wings(1);
            cg=obj.cg;
            fid=fopen(string,'w');
            fprintf(fid,'%s \n',obj.name);
            fprintf(fid,'#Mach \n');
            fprintf(fid,'0.0 \n');
            fprintf(fid,'#IYsym IZsym Zsym \n');
            fprintf(fid,'0 0 0 \n');
            fprintf(fid,'#Sref Cref Bref \n');
            wing_area=0;
            wing_span=0;
            for nw=1:length(obj.wings);
                wing_area=wing_area+obj.wings(nw).s;
                wing_span=wing_span+obj.wings(nw).b;
            end
            fprintf(fid,'%f %f %f \n',wing_area,wing.mean_chord,wing_span);
            fprintf(fid,'#Xref Yref Zref \n');
            fprintf(fid,'%f %f %f \n',cg(1),cg(2),cg(3));
            
            for k=1:length(obj.wings)
                obj.wings(k).build_surface(fid)
            end
            for k=1:length(obj.htails)
                obj.htails(k).build_surface(fid)
            end
            for k=1:length(obj.vtails)
                obj.vtails(k).build_surface(fid)
            end
            fprintf(fid,'#========================================== \n');
            fid=fclose(fid);
        end
        %% BUILD MASS FILE
        function build_mass(obj)
            string=['.\AVL_FILES\',obj.name,'.mass'];
            fid=fopen(string,'wt');
            fid=fclose(fid);
            wing=obj.wings(1);
            cg=obj.cg;
            fid=fopen(string,'w');
            fprintf(fid,'Lunit = 1.0 m\n');
            fprintf(fid,'Munit = 1.0 kg\n');
            fprintf(fid,'Tunit = 1.0 s \n');
            fprintf(fid,'g = %.3f \n',obj.param.g);
            fprintf(fid,'rho = %.3f \n',obj.param.d);
            fprintf(fid,'# mass x y z Ixx Iyy Izz \n');
            fprintf(fid,'%f %f %f %f %f %f %f \n',obj.m,obj.cg(1),obj.cg(2),obj.cg(3),obj.inert(1),obj.inert(2),obj.inert(3));
            fprintf(fid,'#========================================== \n');
            fid=fclose(fid);
        end
        %% RECOVER AVL DATA/Run AVL
        function obj=run_avl(obj)
            filename=['.\AVL_FILES\',obj.name];
            basename=['.\AVL_DATA\',obj.name,'_Data'];
            runAVL(filename,basename,obj);
            st_file=['.\AVL_DATA\',obj.name,'_Data.st'];
            run_file=['.\AVL_DATA\',obj.name,'_Data.sb'];
            eig_file=['.\AVL_DATA\',obj.name,'_Data.eig'];
            obj.st=parseST(st_file);
            obj.run=parseRunCaseHeader(run_file);
            eig_full=importdata(eig_file);
            eig_dat=eig_full.data(:,2:3);
            i=1;
            n=0;
            while n<5
                n=n+1;
                roots(n,:)=eig_dat(i,:);
                if i < length(eig_dat(:,1))
                    if eig_dat(i,1) == eig_dat(i+1,1)
                        i=i+1;
                    end
                end
                i=i+1;
            end
            
            function [zeta,omega]=roots_to_damping(roots)
                real = roots(1);
                imag=roots(2);
                if real < 0
                    zeta=sqrt(1/(1+(imag/real)^2));
                    omega=-1*real/zeta;
                elseif real > 0
                    zeta = -1;
                    omega=0;
                else
                    zeta=0;
                    omega=imag;
                end
            end
            
            [obj.eig.roll(1,1),obj.eig.roll(1,2)]=roots_to_damping(roots(1,:));
            [obj.eig.dutch(1,1),obj.eig.dutch(1,2)]=roots_to_damping(roots(2,:));
            [obj.eig.short(1,1),obj.eig.short(1,2)]=roots_to_damping(roots(3,:));
            [obj.eig.spiral(1,1),obj.eig.spiral(1,2)]=roots_to_damping(roots(4,:));
            [obj.eig.phugoid(1,1),obj.eig.phugoid(1,2)]=roots_to_damping(roots(5,:));
        end
        
    end
end