function bodies = importMultipleBodies(object_files, ...
                                        rotation_hinge_points_CAD, ...
                                        rotation_directions_CAD, ...
                                        temperatures__K, ...
                                        energy_accommodation_coefficients, ...
                                        DCM_B_from_CAD, ...
                                        CoM_CAD)
%% importMultipleBodies - Imports multiple bodies from .obj files and prepares them for the VLEO aerodynamics simulation
%
%  bodies = importMultipleBodies(object_files, ...
%                                rotation_hinge_points_CAD, ...
%                                rotation_directions_CAD, ...
%                                temperatures__K, ...
%                                energy_accommodation_coefficients, ...
%                                DCM_B_from_CAD, ...
%                                CoM_CAD)
%
%  This function imports multiple bodies from .obj files and prepares them
%  for the VLEO aerodynamics simulation. The function returns a cell array
%  of structures containing the vertices, surface centroids, surface normals,
%  rotation direction, rotation hinge point, surface temperatures,
%  surface energy accommodation coefficients, and surface areas of the bodies.
%
%  Inputs:
%   object_files: 1xN array of strings of the paths to the .obj or .m (gmsh) files
%   rotation_hinge_points_CAD: 3xN array of the rotation hinge points of the bodies in the CAD frame
%   rotation_directions_CAD: 3xN array of the rotation directions of the bodies in the CAD frame
%   temperatures__K: 1xN cell array of the surface temperatures of the bodies
%   energy_accommodation_coefficients: 1xN cell array of the surface energy accommodation coefficients of the bodies
%   DCM_B_from_CAD: 3x3 array of the direction cosine matrix from the CAD frame to the body frame
%   CoM_CAD: 3x1 array of the center of mass of the bodies in the CAD frame
%
%  Outputs:
%   bodies: 1xN cell array of structures containing the vertices, surface centroids, surface normals,
%           rotation direction, rotation hinge point, surface temperatures,
%           surface energy accommodation coefficients, and surface areas of the bodies
%

arguments
    object_files (1,:) string {mustBeFile}
    rotation_hinge_points_CAD (3,:) {mustBeNumeric, mustBeReal}
    rotation_directions_CAD (3,:) {mustBeNumeric, mustBeReal}
    temperatures__K (1,:) cell {}
    energy_accommodation_coefficients (1,:) cell {}
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
% Check if first file is .gmsh and validate accordingly
[~, ~, first_ext] = fileparts(object_files(1));
if strcmpi(first_ext, '.m')
    if length(object_files) > 1
        error('When using .gmsh files, only one file is allowed.');
    end
    body_data = import_gmsh(object_files(1));
end

% Prepare structure of objects
num_bodies = size(rotation_hinge_points_CAD,2);
bodies = cell(1, num_bodies);

for i = 1:num_bodies
    switch lower(first_ext)
        case '.obj'
            %% Extract data from obj file
            current_object_file = object_files(i);

            % Vertices
            vertices_CAD = getVerticesFromObj(current_object_file);
        case '.m'
            vertices_CAD = body_data(i).vertices_CAD; % These are in CAD frame
        otherwise
            error('Unsupported file format: %s. Only .obj and .gmsh files are supported.', first_ext);
    end
    % Transform vertices to body frame
    vertices_CAD_list = reshape(vertices_CAD, 3, []);
    vertices_CAD_list = vertices_CAD_list - CoM_CAD;
    vertices_B_list = DCM_B_from_CAD * vertices_CAD_list;
    vertices_B = reshape(vertices_B_list, 3, 3, []);
    bodies{i}.vertices_B = vertices_B;

    % Calculate surface centroids
    bodies{i}.centroids_B = squeeze(mean(vertices_B, 2));

    % Calculate surface normals
    AB_vectors = squeeze(vertices_B(:,2,:) - vertices_B(:,1,:));
    AC_vectors = squeeze(vertices_B(:,3,:) - vertices_B(:,1,:));
    cross_products = cross(AB_vectors, AC_vectors);
    cross_product_norms = vecnorm(cross_products);
    bodies{i}.normals_B = cross_products ./ cross_product_norms;

    % Calculate surface areas
    bodies{i}.areas = 1/2 * cross_product_norms;

    %% Write remaining data to bodies structure

    % Rotation hinge point
    bodies{i}.rotation_hinge_point_B = DCM_B_from_CAD * (rotation_hinge_points_CAD(:, i) - CoM_CAD);
    % Rotation direction
    bodies{i}.rotation_direction_B = DCM_B_from_CAD * rotation_directions_CAD(:, i);

    % Surface temperatures and energy accommodation coefficients
    current_num_faces = size(vertices_B, 3);

    current_temperatures__K = temperatures__K{i};
    bodies{i}.temperatures__K = nan(1, current_num_faces);
    bodies{i}.temperatures__K(:) = current_temperatures__K; % if scalar, it will be expanded to the correct size

    current_energy_accommodation_coefficients = energy_accommodation_coefficients{i};
    bodies{i}.energy_accommodation_coefficients = nan(1, current_num_faces);
    bodies{i}.energy_accommodation_coefficients(:) = current_energy_accommodation_coefficients; % if scalar, it will be expanded to the correct size
end
end