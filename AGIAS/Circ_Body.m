classdef Circ_Body
    properties
        l %length of body
        s_array %shape array
        coord %translates body [X,Y,Z]
    end
    
    methods
        function obj=Circ_Body(L,S_array,Coord)
            obj.l=L;
            obj.s_array=S_array;
            obj.coord=Coord;
        end
        function plot_body(obj)
            hold on
            n=length(obj.s_array);
            theta=linspace(0,2*pi,180);
            x=linspace(0,obj.l,n);
            
            for i=1:n
                Y=obj.s_array(i)*cos(theta)+obj.coord(2);
                Z=obj.s_array(i)*sin(theta)+obj.coord(3);
                X=Z.*0+x(i)+obj.coord(1);
                plot3(X,Y,Z,'k')
            end
            axis equal
        end
    end
end
