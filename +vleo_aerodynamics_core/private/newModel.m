function [aeroForce__N, aeroTorque__Nm] = newModel(areas__m2,...
    normals,...
    centroids__m,...
    v_rels__m_per_s,...
    deltas__rad,...
    density__kg_per_m3,...
    aerodynamic_coefficents)
    %% newModel - computes aerodynamic forces based on the new IRS Model.
    % Inputs:
    %   areas__m2: 1xN array of the areas of N triangles
    %   normals: 3xN array of surface normals of N triangles
    %   centroids__m: 3xN array of surface centroids of N triangles
    %   v_rels__m_per_s: 3xN array of relative velocities of N triangles
    %   deltas__rad: 1xN array of angles between the flow direction and the normals of N triangles
    %   density__kg_per_m3: scalar value of the incomming streams density
    %   aerodynamic_coefficent_functions: struct with one field for C_l and
    %                                     C_d in dependece of Angle of
    %                                     Attack:
    %                                     1. field name: curve_c_l
    %                                     2. field name: curve_c_d
    %   aerodynamic_coefficent_functions: a 'griddedInterpolant' object containing the lookup table data for
    %             the 2 aerodynamic coefficients,negative angles of attack resemble wake faces:
    %             -  C_l
    %             -  C_d
    % Outputs:
    %   aeroForce__N: 3x1 array of the aerodynamic force acting on the body in the same coordinate
    %                 system as the inputs normals and centroids
    %   aeroTorque__Nm: 3x1 array of the aerodynamic torque acting on the body in the same
    %                   coordinate system as the inputs normals and centroids and with respect to its origin 
    %
    arguments
        areas__m2 (1,:) {mustBeNumeric, mustBeReal, mustBePositive};
        normals (3,:) {mustBeNumeric, mustBeReal};
        centroids__m (3,:) {mustBeNumeric, mustBeReal};
        v_rels__m_per_s (3,:) {mustBeNumeric, mustBeReal};
        deltas__rad (1,:) {mustBeNumeric, mustBeReal};
        density__kg_per_m3 (1,1) {mustBeNumeric, mustBeReal, mustBePositive};
        aerodynamic_coefficents ;
    end
    %% assertions
    if isstruct(aerodynamic_coefficents) 
        % required_fields = {'curve_c_l','curve_c_d'};
        % for i = 1:numel(required_fields)
        %     if ~isfield(aerodynamic_coefficents,required_fields{i})
        %         error('Missing required field: %s', required_fields{i});
        %     end
        % end
    elseif isa(aerodynamic_coefficents,'griddedInterpolant')
        % Replace assert with if-statement to avoid compile-time evaluation
        if size(aerodynamic_coefficents.Values, 2) ~= 2
            error('LUT_data must return a matrix with 2 columns for C_l, C_d');
        end
    else
        error('aerodynamic_coefficent_functions must be a struct with fields curve_c_l and curve_c_d or a griddedInterpolant object.');
    end
    
    %% Abbreviations
    v_rels = v_rels__m_per_s;
    V = vecnorm(v_rels);
    rho = density__kg_per_m3;
    v_hat = v_rels./V;

    %%LUT
    AOA__deg = 90-deltas__rad*180/pi;
    if isstruct(aerodynamic_coefficents)
        % Use the struct with fields curve_c_l and curve_c_d
        C_l = ppval(aerodynamic_coefficents.pp_Cl,AOA__deg);
        C_d = ppval(aerodynamic_coefficents.pp_Cd,AOA__deg);
    else
        c = aerodynamic_coefficents(AOA__deg);
        C_l = c(:,1)';
        C_d = c(:,2)';
    end

    %darg
    F_d_mag = 0.5*rho*V.^2.*areas__m2.*C_d;
    F_l_mag = 0.5*rho*V.^2.*areas__m2.*C_l;
    F_d = F_d_mag.*v_hat;

    %lift
    lift_dir = -cross(cross(v_hat,normals),v_hat);
    lift_dir_norm = vecnorm(lift_dir);
    lift_dir_norm(lift_dir_norm == 0) = 1; % avoid division by zero
    F_l = F_l_mag.*lift_dir./lift_dir_norm;

    %resultant force
    aeroForce__N = F_l + F_d;
    aeroTorque__Nm = cross(centroids__m,aeroForce__N,1);

end