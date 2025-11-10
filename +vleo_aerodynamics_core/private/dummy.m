function [aeroForce__N, aeroTorque__Nm] = dummy(areas__m2,...
                                                v_rels__m_per_s,...
                                                centroids__m,...
                                                deltas__rad)

%DUMMY model for test purposes that always outputs 1N of Drag
%   For testing purposes this model will use 1N/m^2 of drag for flow
%   facing triangles, not flow facing faces will have no drag.
%   and calculate the torques based on that drag
%
% Inputs:
%   areas__m2: 1xN array of the areas of N triangles
%   centroids__m: 3xN array of surface centroids of N triangles
%   v_rels__m_per_s: 3xN array of relative velocities of N triangles
%   deltas__rad: 1xN array of angles between the flow direction and the normals of N triangles
%
% Outputs:
%   aeroForce__N: 3x1 array of the aerodynamic force acting on the body in the same coordinate
%                 system as the inputs normals and centroids
%   aeroTorque__Nm: 3x1 array of the aerodynamic torque acting on the body in the same
%                   coordinate system as the inputs normals and centroids and with respect to its origin 
%
%% Abbreviations
v_rels = v_rels__m_per_s;
V = vecnorm(v_rels);
v_hat = (v_rels./V).*areas__m2;

%% detect backward facing faces
backward_facing = deltas__rad >= pi/2;
aeroForce__N= v_hat;
aeroForce__N(:,backward_facing) = 0;
aeroTorque__Nm = cross(centroids__m,aeroForce__N);
end

