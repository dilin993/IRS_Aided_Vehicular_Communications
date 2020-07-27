classdef node
    
    properties
        id;
        Nx;
        Ny;
        Nz;
        N;
        delta;
        pos;
        arrival_phi;
        departure_phi;
        arrival_theta;
        departure_theta;
        gain;
    end
    
    methods
        function obj = node(id, Nx, Ny, Nz, delta, pos, gain)
            obj.id = id;
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.Nz = Nz;
            obj.N = Nx*Ny*Nz;
            obj.delta = delta;
            obj.pos = pos;
            obj.arrival_phi = -pi/2 + pi * rand;
            obj.arrival_theta = -pi/4 + (pi/2) * rand;
            obj.departure_phi = -pi/2 + pi * rand;
            obj.departure_theta = -pi/4 + (pi/2) * rand;
            obj.gain = gain;
        end
        
        function obj = change_angles(obj)
            obj.arrival_phi = -pi/2 + pi * rand;
            obj.arrival_theta = -pi/4 + (pi/2) * rand;
            obj.departure_phi = -pi/2 + pi * rand;
            obj.departure_theta = -pi/4 + (pi/2) * rand;
        end
    end
end

