function showBodies(bodies, opts)
% showBodies - Plot the bodies and their surface centroids and normals rotated by the given angles
%
%  showBodies(bodies, opts)
%
%  This function plots the bodies and their surface normals rotated by the given angles.
%  The bodies are plotted in a 3D figure with the surface centroids and normals. In addition,
%  scalar and vectorial values can be given to be plotted on the surfaces. Scalar values are
%  plotted as surface colors and vectorial values are plotted as vectors.
%
%  Inputs:
%   bodies: 1xN cell array of structures containing the vertices, surface centroids, surface normals,
%           rotation direction, rotation hinge point
%   opts: Name-value arguments structure with optional fields
%

arguments
    bodies (1,:) cell
    opts.bodies_rotation_angles (1,:) double {mustBeReal} = 0
    opts.face_alpha (1,1) {mustBeNonnegative, mustBeLessThanOrEqual(opts.face_alpha,1)} = 0.75
    opts.scale_normals (1,1) double {mustBeNonnegative} = 1
    opts.scalar_values cell = {}
    opts.vectorial_values cell = {}
    opts.scale_vectorial_values (1,1) double {mustBeNonnegative} = 1
end

num_bodies = length(bodies);

% Handle all parameters with consistent if-else logic
if length(opts.bodies_rotation_angles) == 1 && opts.bodies_rotation_angles == 0
    rotation_angles = zeros(1, num_bodies);
else
    rotation_angles = opts.bodies_rotation_angles;
end

if isempty(opts.scalar_values)
    scalar_values = cell(1, num_bodies);
else
    scalar_values = opts.scalar_values;
end

if isempty(opts.vectorial_values)
    vectorial_values = cell(1, num_bodies);
else
    vectorial_values = opts.vectorial_values;
end

figure;
hold on;
axis equal;
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
if ~isempty(scalar_values) && any(~cellfun(@isempty, scalar_values))
    colorbar;
end

for i = 1:num_bodies
    current_body = bodies{i};
    current_rotation_angle = rotation_angles(i);
    current_scalar_values = scalar_values{i};
    current_vectorial_values = vectorial_values{i};

    % Rotate the body
    [vertices_B, centroids_B, normals_B] = rotateBody(current_body.vertices_B, ...
                                                     current_body.centroids_B, ...
                                                     current_body.normals_B, ...
                                                     current_rotation_angle, ...
                                                     current_body.rotation_direction_B, ...
                                                     current_body.rotation_hinge_point_B);

    X = squeeze(vertices_B(1, :, :));
    Y = squeeze(vertices_B(2, :, :));
    Z = squeeze(vertices_B(3, :, :));

    num_faces = size(X, 2);
    
    % Handle scalar values for coloring
    if isempty(current_scalar_values)
        C = 'green';
    else
        if length(current_scalar_values) == num_faces
            C = current_scalar_values;
        elseif length(current_scalar_values) == 1
            C = repmat(current_scalar_values, 1, num_faces);
        else
            warning('Scalar values size mismatch for body %d. Using default color.', i);
            C = 'green';
        end
    end
    
    % Plot the body
    patch(X, Y, Z, C, 'FaceAlpha', opts.face_alpha, 'DisplayName', sprintf('Body %d', i));

    % Plot hinge point, rotation direction, centroids, and normals - always
    scatter3(current_body.rotation_hinge_point_B(1), current_body.rotation_hinge_point_B(2), current_body.rotation_hinge_point_B(3), 'filled', 'black');
    quiver3(current_body.rotation_hinge_point_B(1), current_body.rotation_hinge_point_B(2), current_body.rotation_hinge_point_B(3), ...
            current_body.rotation_direction_B(1), current_body.rotation_direction_B(2), current_body.rotation_direction_B(3), 'Color', 'black');
    scatter3(centroids_B(1, :), centroids_B(2, :), centroids_B(3, :), [], 'red', 'filled');
    quiver3(centroids_B(1, :), centroids_B(2, :), centroids_B(3, :), ...
            normals_B(1, :), normals_B(2, :), normals_B(3, :), opts.scale_normals, 'Color', 'red');

    % Plot vectorial values if provided
    if ~isempty(current_vectorial_values)
        if size(current_vectorial_values, 2) == size(centroids_B, 2)
            quiver3(centroids_B(1, :), centroids_B(2, :), centroids_B(3, :), ...
                current_vectorial_values(1, :), current_vectorial_values(2, :), current_vectorial_values(3, :), ...
                opts.scale_vectorial_values, 'Color', 'b');
        else
            warning('Vectorial values size mismatch for body %d. Skipping.', i);
        end
    end
end

legend('Face', 'Rotation Hinge Point', 'Rotation Direction', 'Face Centroid', 'Face Normal')

end

% % showBodies - Plot the bodies and their surface centroids and normals rotated by the given angles
% %
% %  showBodies(bodies, bodies_rotation_angles__rad, face_alpha, scale_normals, scalar_values, vectorial_values, scale_vectorial_values)
% %
% %  This function plots the bodies and their surface normals rotated by the given angles.
% %  The bodies are plotted in a 3D figure with the surface centroids and normals. In addition,
% %  scalar and vectorial values can be given to be plotted on the surfaces. Scalar values are
% %  plotted as surface colors and vectorial values are plotted as vectors.
% %
% %  Inputs:
% %   bodies: 1xN cell array of structures containing the vertices, surface centroids, surface normals,
% %           rotation direction, rotation hinge point
% %   bodies_rotation_angles__rad: 1xN array of the rotation angles of the bodies
% %   face_alpha: scalar, the transparency of the faces
% %   scale_normals: scalar, the scale factor for the normals
% %   scalar_values: 1xN cell array of scalar values to be plotted on the surfaces
% %   vectorial_values: 1xN cell array of vectorial values to be plotted on the surfaces
% %   scale_vectorial_values: scalar, the scale factor for the vectorial values
% %
% 
% arguments
%     bodies (1,:) cell
%     bodies_rotation_angles__rad (1,:) double {mustBeReal}
%     face_alpha (1,1) {mustBeNonnegative, mustBeLessThanOrEqual(face_alpha,1)} = 0.75
%     scale_normals double {mustBeNonnegative} = 1
%     scalar_values cell = cell(1, length(bodies))
%     vectorial_values cell = cell(1, length(bodies))
%     scale_vectorial_values double {mustBeNonnegative} = 1
% end
% 
% figure;
% hold on;
% axis equal;
% grid on;
% xlabel('x');
% ylabel('y');
% zlabel('z');
% if ~isempty(scalar_values)
%     colorbar;
% end
% 
% num_bodies = length(bodies);
% for i = 1:num_bodies
%     current_body = bodies{i};
%     current_rotation_angle__rad = bodies_rotation_angles__rad(i);
%     current_scalar_values = scalar_values{i};
%     current_vectorial_values = vectorial_values{i};
% 
%     % Rotate the body
%     [vertices_B, centroids_B, normals_B] = rotateBody(current_body.vertices_B, ...
%                                                              current_body.centroids_B, ...
%                                                              current_body.normals_B, ...
%                                                              current_rotation_angle__rad, ...
%                                                              current_body.rotation_direction_B, ...
%                                                              current_body.rotation_hinge_point_B);
% 
%     X = squeeze(vertices_B(1, :, :));
%     Y = squeeze(vertices_B(2, :, :));
%     Z = squeeze(vertices_B(3, :, :));
% 
%     % Plot the body
%     if isempty(current_scalar_values)
%         C = 'green';
%     else
%         C = current_scalar_values;
%     end
%     patch(X, Y, Z, C, 'FaceAlpha', face_alpha, 'DisplayName', sprintf('Body %d', i));
% 
%     % Plot the hinge point
%     scatter3(current_body.rotation_hinge_point_B(1), ...
%              current_body.rotation_hinge_point_B(2), ...
%              current_body.rotation_hinge_point_B(3), 'filled', 'black');
% 
%     % Plot the rotation direction
%     quiver3(current_body.rotation_hinge_point_B(1), ...
%             current_body.rotation_hinge_point_B(2), ...
%             current_body.rotation_hinge_point_B(3), ...
%             current_body.rotation_direction_B(1), ...
%             current_body.rotation_direction_B(2), ...
%             current_body.rotation_direction_B(3), 'Color', 'black');
% 
%     % Plot the centroids
%     scatter3(centroids_B(1, :), centroids_B(2, :), centroids_B(3, :), [], 'red', 'filled');
% 
%     % Plot the normals
%     quiver3(centroids_B(1, :), centroids_B(2, :), centroids_B(3, :), ...
%             normals_B(1, :), normals_B(2, :), normals_B(3, :), scale_normals,'Color', 'red');
% 
%     % Plot the vectorial values
%     if ~isempty(current_vectorial_values)
%         quiver3(centroids_B(1, :), centroids_B(2, :), centroids_B(3, :), ...
%             current_vectorial_values(1, :), current_vectorial_values(2, :), current_vectorial_values(3, :), ...
%             scale_vectorial_values, 'Color', 'b');
% 
%     end
% 
% end
% 
% legend('Face', 'Rotation Hinge Point', 'Rotation Direction', 'Face Centroid', 'Face Normal')
% 
% end