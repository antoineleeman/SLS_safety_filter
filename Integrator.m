
classdef Integrator

    properties
        nx = 2; % number of state
        nu = 1; % number of input
        ni = 6;
        nw = 2;
        
        x0 = [0;0];
        xf = [0;0];
        
        A = [1,1;0,1];
        B = [0.5;1];
        
        F_x;
        b_x; 
        
        F_u;
        b_u; 
        
        Bw;
        
        N = 10;
        dt;
        
        Q_cost;
        R_cost;
        
        x_max;
        u_max;
        
    end
    
    methods
        function obj = Integrator()
            obj.Bw = 0.3*eye(obj.nx);
            obj.Q_cost = eye(obj.nx);
            obj.R_cost = 100*eye(obj.nu);
            obj.x_max = 5;
            obj.u_max = 3;
            C = [0,0,-1;0,0,1;1,0,0;-1,0,0;0,1,0;0,-1,0];
            c = [obj.u_max ;obj.u_max ;obj.x_max;obj.x_max;obj.x_max;obj.x_max];
            obj.F_x = C(3:end,1:obj.nx);
            obj.b_x = c(3:end);
            
            obj.F_u = C(1:2,obj.nx+1:end);
            obj.b_u = c(1:2);
        end
        
        
    end
end

