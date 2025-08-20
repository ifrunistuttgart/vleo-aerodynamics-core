function body_data = import_obj(obj_file, DCM_B_from_CAD, CoM_CAD)
%% import_obj - Import a single .obj file and transform to body frame
%
%  body_data = import_obj(obj_file, DCM_B_from_CAD, CoM_CAD)
%
%  This function imports a single .obj file and transforms the vertices
%  from CAD frame to body frame, then calculates centroids, normals, and areas.
%
%  Inputs:
%   obj_file: string path to the .obj file
%   DCM_B_from_CAD: 3x3 direction cosine matrix from CAD frame to body frame
%   CoM_CAD: 3x1 center of mass in CAD frame
%
%  Outputs:
%   body_data: structure containing:
%     - vertices_B: 3x3xN array of triangle vertices in body frame
%     - centroids_B: 3xN array of triangle centroids in body frame
%     - normals_B: 3xN array of triangle normal vectors in body frame
%     - areas: 1xN array of triangle areas

    % Extract vertices from obj file
    vertices_CAD = getVerticesFromObj(obj_file);
    
    % Transform vertices to body frame
    vertices_CAD_list = reshape(vertices_CAD, 3, []);
    vertices_CAD_list = vertices_CAD_list - CoM_CAD;
    vertices_B_list = DCM_B_from_CAD * vertices_CAD_list;
    vertices_B = reshape(vertices_B_list, 3, 3, []);
    
    % Calculate surface centroids
    centroids_B = squeeze(mean(vertices_B, 2));
    
    % Calculate surface normals
    AB_vectors = squeeze(vertices_B(:,2,:) - vertices_B(:,1,:));
    AC_vectors = squeeze(vertices_B(:,3,:) - vertices_B(:,1,:));
    cross_products = cross(AB_vectors, AC_vectors);
    cross_product_norms = vecnorm(cross_products);
    normals_B = cross_products ./ cross_product_norms;
    
    % Calculate surface areas
    areas = 1/2 * cross_product_norms;   %areas instead of areas_B
    
    % Create output structure
    body_data = struct('vertices_B', vertices_B, ...
                       'centroids_B', centroids_B, ...
                       'normals_B', normals_B, ...
                       'areas', areas);
end