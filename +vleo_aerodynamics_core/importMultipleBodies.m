function bodies = importMultipleBodies(filepaths, ...
                                        rotation_hinge_points_CAD, ...
                                        rotation_directions_CAD, ...
                                        temperatures__K, ...
                                        energy_accommodation_coefficients, ...
                                        DCM_B_from_CAD, ...
                                        CoM_CAD)
%% importMultipleBodies - Imports multiple bodies from .obj files or a single .m file and prepares them for the VLEO aerodynamics simulation
% %
% %  bodies = importMultipleBodies(filepaths, ...
% %                                rotation_hinge_points_CAD, ...
% %                                rotation_directions_CAD, ...
% %                                temperatures__K, ...
% %                                energy_accommodation_coefficients, ...
% %                                DCM_B_from_CAD, ...
% %                                CoM_CAD)
% %
% %  This function imports multiple bodies from .obj files and prepares them
% %  for the VLEO aerodynamics simulation. The function returns a cell array
% %  of structures containing the vertices, surface centroids, surface normals,
% %  rotation direction, rotation hinge point, surface temperatures,
% %  surface energy accommodation coefficients, and surface areas of the bodies.
% %
% %  Inputs:
% %   filepaths: 1xN array of strings of the paths to the .obj or .m files
% %   rotation_hinge_points_CAD: 3xN array of the rotation hinge points of the bodies in the CAD frame
% %   rotation_directions_CAD: 3xN array of the rotation directions of the bodies in the CAD frame
% %   temperatures__K: 1xN cell array of the surface temperatures of the bodies
% %   energy_accommodation_coefficients: 1xN cell array of the surface energy accommodation coefficients of the bodies
% %   DCM_B_from_CAD: 3x3 array of the direction cosine matrix from the CAD frame to the body frame
% %   CoM_CAD: 3x1 array of the center of mass of the bodies in the CAD frame
% %
% %  Outputs:
% %   bodies: 1xN cell array of structures containing the vertices, surface centroids, surface normals,
% %           rotation direction, rotation hinge point, surface temperatures,
% %           surface energy accommodation coefficients, and surface areas of the bodies
% %
arguments
    filepaths (1,:) string {mustBeFile}
    rotation_hinge_points_CAD (3,:) {mustBeNumeric, mustBeReal}
    rotation_directions_CAD (3,:) {mustBeNumeric, mustBeReal}
    temperatures__K (1,:) cell
    energy_accommodation_coefficients (1,:) cell
    DCM_B_from_CAD (3,3) {mustBeNumeric, mustBeReal}
    CoM_CAD (3,1) {mustBeNumeric, mustBeReal}
end

% Check if DCM is orthogonal
D = det(DCM_B_from_CAD);
scaled_DCM = DCM_B_from_CAD / nthroot(D,3);
if norm(scaled_DCM * scaled_DCM' - eye(3)) > 1e-6
    error('DCM is not orthogonal.');
end
if D ~= 1
    warning('Determinant of DCM is not equal to 1. The DCM will be scaled accordingly.');
    DCM_B_from_CAD = scaled_DCM;
    disp('New DCM_B_from_CAD:');
    disp(DCM_B_from_CAD);
end

% Determine file type and import accordingly
[~, ~, exts] = fileparts(filepaths);
unique_exts = unique(exts);

if length(unique_exts) == 1 && strcmpi(unique_exts{1}, '.obj')
    fprintf('Processing multiple .obj files...\n');
    bodies = cell(1, length(filepaths));
    
    for i = 1:length(filepaths)
        % Import each .obj file
        body_data = import_obj(filepaths(i), DCM_B_from_CAD, CoM_CAD);
        bodies{i} = body_data;
    end
    
elseif length(filepaths) == 1 && strcmpi(exts{1}, '.m')
    fprintf('Processing single .m file with multiple bodies...\n');
    bodies = import_gmsh(filepaths(1));
    
    % Transform all bodies to body frame
    for i = 1:length(bodies)
        vertices_CAD = bodies{i}.vertices_B;
        vertices_CAD_list = reshape(vertices_CAD, 3, []);
        vertices_CAD_list = vertices_CAD_list - CoM_CAD;
        vertices_B_list = DCM_B_from_CAD * vertices_CAD_list;
        bodies{i}.vertices_B = reshape(vertices_B_list, 3, 3, []);
        
        % Recalculate centroids in body frame
        bodies{i}.centroids_B = squeeze(mean(bodies{i}.vertices_B, 2));
        
        % Recalculate normals in body frame
        AB_vectors = squeeze(bodies{i}.vertices_B(:,2,:) - bodies{i}.vertices_B(:,1,:));
        AC_vectors = squeeze(bodies{i}.vertices_B(:,3,:) - bodies{i}.vertices_B(:,1,:));
        cross_products = cross(AB_vectors, AC_vectors);
        cross_product_norms = vecnorm(cross_products);
        bodies{i}.normals_B = cross_products ./ cross_product_norms;
        
    end
    
else
    error('Unsupported file configuration. Either provide multiple .obj files or a single .m file.');
end

% Get the actual number of bodies
num_bodies = length(bodies);

% Add remaining data to bodies structure
for i = 1:num_bodies
    % Rotation hinge point
    bodies{i}.rotation_hinge_point_B = DCM_B_from_CAD * (rotation_hinge_points_CAD(:, i) - CoM_CAD);
    
    % Rotation direction
    bodies{i}.rotation_direction_B = DCM_B_from_CAD * rotation_directions_CAD(:, i);
    
    % Surface temperatures and energy accommodation coefficients
    current_num_faces = size(bodies{i}.vertices_B, 3);
    
    current_temperatures__K = temperatures__K{i};
    bodies{i}.temperatures__K = nan(1, current_num_faces);
    bodies{i}.temperatures__K(:) = current_temperatures__K;
    
    current_energy_accommodation_coefficients = energy_accommodation_coefficients{i};
    bodies{i}.energy_accommodation_coefficients = nan(1, current_num_faces);
    bodies{i}.energy_accommodation_coefficients(:) = current_energy_accommodation_coefficients; % if scalar, it will be expanded to the correct size
end
end
