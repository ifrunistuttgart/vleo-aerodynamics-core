function [aerodynamic_force_B__N, ...
            aerodynamic_torque_B_B__Nm] = ...
                            simplifiedVleoAerodynamics(attitude_quaternion_BI, ...
                                                        rotational_velocity_BI_B__rad_per_s, ...
                                                        velocity_I_I__m_per_s, ...
                                                        wind_velocity_I_I__m_per_s, ...
                                                        density__kg_per_m3, ...
                                                        temperature__K, ... 
                                                        particles_mass__kg, ...
                                                        bodies, ...                                                       
                                                        bodies_rotation_angles__rad, ...
                                                        temperature_ratio_method)
%% simplifiedVleoAerodynamics - Calculate the aerodynamic force and torque acting on a satellite in VLEO
%
%   [aeroForce__N, aeroTorque__Nm] = simplifiedVleoAerodynamics(attitude_quaternion_BI, ...
%                                                               rotational_velocity_BI_B__rad_per_s, ...
%                                                               velocity_I_I__m_per_s, ...
%                                                               wind_velocity_I_I__m_per_s, ...
%                                                               density__kg_per_m3, ...
%                                                               temperature__K, ...
%                                                               particles_mass__kg, ...
%                                                               bodies, ...
%                                                               bodies_rotation_angles__rad, ...
%                                                               temperature_ratio_method)
%
%   This function calculates the aerodynamic force and torque acting on a satellite in VLEO.
%
%   This function expects the space_math_utilities namespace to be available on the MATLAB path.
%
%   Inputs:
%    attitude_quaternion_BI: 4x1 array of the attitude quaternion from the body frame to the inertial frame
%    rotational_velocity_BI_B__rad_per_s: 3x1 array of the rotational velocity of the satellite with respect
%                                         to the inertial frame expressed in the body frame
%    velocity_I_I__m_per_s: 3x1 array of the velocity of the satellite with respect to the inertial frame
%                           expressed in the inertial frame
%    wind_velocity_I_I__m_per_s: 3x1 array of the velocity of the wind with respect to the inertial frame
%                                expressed in the inertial frame
%    density__kg_per_m3: Scalar value of the density of the gas
%    temperature__K: Scalar value of the temperature of the gas
%    particles_mass__kg: Scalar value of the mass of the particles
%    bodies: 1xN cell array of structures containing the vertices, surface centroids, surface normals,
%            rotation direction, rotation hinge point, surface temperatures,
%            surface energy accommodation coefficients, and surface areas of the bodies
%    bodies_rotation_angles__rad: 1xN array of the rotation angles of the bodies
%    temperature_ratio_method: Scalar value of the method to calculate the temperature ratio
%                              1: Exact term
%                              2: Hyperthermal approximation 1
%                              3: Hyperthermal approximation 2
%
%  Outputs:
%   aerodynamic_force_B__N: 3x1 array of the aerodynamic force acting on the satellite expressed in the body frame
%   aerodynamic_torque_B_B__Nm: 3x1 array of the aerodynamic torque acting on the satellite with respect to the 
%                               center of mass (origin of body frame) expressed in the body frame
%

import space_math_utilities.cpm

%% Abbreviations
q_BI = quaternion(attitude_quaternion_BI');
omega = rotational_velocity_BI_B__rad_per_s;
v_rel_I = wind_velocity_I_I__m_per_s - velocity_I_I__m_per_s;

%% Extract data from bodies structure
% Get total number of faces of all bodies
num_bodies = length(bodies);
bodies_num_faces = zeros(1, num_bodies);
for i = 1:num_bodies
    current_body = bodies{i};
    bodies_num_faces(i) = size(current_body.vertices_B,3);
end
total_num_faces = sum(bodies_num_faces);

% Prepare variables
% Surface temperatures
surface_temperatures__K = nan(1, total_num_faces);
% Energy accommodation coefficients
energy_accommodation_coefficients = surface_temperatures__K;
% Vertex coordinates
vertices_B = nan(3, 3, total_num_faces);
% Surface normals
normals_B = nan(3, total_num_faces);
% Surface centroids
centroids_B = normals_B;
% Surface Areas
areas = surface_temperatures__K;

last_face_idx = 0;
for i = 1:num_bodies
    current_body = bodies{i};

    current_indices = last_face_idx + (1:bodies_num_faces(i));

    %% Extract nondirectional data from current body
    surface_temperatures__K(current_indices) = current_body.temperatures__K;
    energy_accommodation_coefficients(current_indices) = current_body.energy_accommodation_coefficients;
    areas(current_indices) = current_body.areas;
    
    %% Rotate directional data according to bodies_rotation_angles__rad  
    current_angle__rad = bodies_rotation_angles__rad(i);
    current_rotation_direction_B = current_body.rotation_direction_B;
    current_rotation_hinge_point_B = current_body.rotation_hinge_point_B;

    [vertices_B(:,:,current_indices), ...
     centroids_B(:,current_indices), ...
     normals_B(:,current_indices)] = rotateBody(current_body.vertices_B, ...
                                                current_body.centroids_B, ...
                                                current_body.normals_B, ...
                                                current_angle__rad, ...
                                                current_rotation_direction_B, ...
                                                current_rotation_hinge_point_B);

    %% Update last face index for the next iteration
    last_face_idx = current_indices(end);
end

%% Determine shadowed faces
%  To calculate the aerodynamics, only the faces that are not shadowed by other faces
%  are considered. The shadowing is only determined by approximation. Effects due
%  the rotational velocity of the satellite or due to the thermal velocity of the
%  residual atmosphere are ignored. 

% Transform wind velocity into body frame
v_rel_B = rotateframe(q_BI, v_rel_I')';
v_rel_dir_B = v_rel_B ./ norm(v_rel_B);

ind_not_shadowed = ~determineShadowedTriangles(vertices_B, centroids_B, normals_B, v_rel_dir_B);

%% Calculate forces and torques
% Determine individual relative velocity of each face by adding the term due to rotation
% of satellite to the relative velocity of the satellite's center of mass
v_indiv_B = v_rel_B - cpm(omega) * centroids_B(:, ind_not_shadowed);

% Direction of individual relative velocity
v_indiv_norm = vecnorm(v_indiv_B);
v_indiv_dir_B = v_indiv_B ./ v_indiv_norm;

% Individual angles between flow and normals
deltas = real(acos(dot(-v_indiv_dir_B, normals_B(:,ind_not_shadowed))));

% Forces and Torques
[aerodynamic_force_B__N, aerodynamic_torque_B_B__Nm] = calcAeroForceAndTorque(areas(ind_not_shadowed), ...
                                    normals_B(:, ind_not_shadowed), ...
                                    centroids_B(:, ind_not_shadowed), ...
                                    v_indiv_B, ...
                                    deltas, ...
                                    density__kg_per_m3, ...
                                    temperature__K, ...
                                    surface_temperatures__K(ind_not_shadowed), ...
                                    energy_accommodation_coefficients(ind_not_shadowed), ...
                                    particles_mass__kg, ...
                                    temperature_ratio_method);

end