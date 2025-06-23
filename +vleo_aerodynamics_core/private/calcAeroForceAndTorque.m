function [aeroForce__N, aeroTorque__Nm] = calcAeroForceAndTorque(areas__m2, ...
                                                                        normals, ...
                                                                        centroids__m, ...
                                                                        v_rels__m_per_s, ...
                                                                        deltas__rad, ...
                                                                        density__kg_per_m3, ...
                                                                        gas_temperature__K, ...
                                                                        surface_temperatures__K, ...
                                                                        energy_accommodation_coefficients, ...
                                                                        particles_mass__kg, ...
                                                                        temperature_ratio_method)
% calcAeroForceAndTorque - Calculate the aerodynamic force and torque acting on a body
%
%   [aeroForce__N, aeroTorque__Nm] = calcAeroForceAndTorque(areas__m2, ...
%                                                           normals, ...
%                                                           centroids__m, ...
%                                                           v_rels__m_per_s, ...
%                                                           deltas__rad, ...
%                                                           density__kg_per_m3, ...
%                                                           gas_temperature__K, ...
%                                                           surface_temperatures__K, ...
%                                                           energy_accommodation_coefficients, ...
%                                                           particles_mass__kg, ...
%                                                           temperature_ratio_method)
%
%   This function calculates the aerodynamic force and torque acting on a body
%   due to the impact of particles on its surface.
%
%   Inputs:
%   areas__m2: 1xN array of the areas of N triangles
%   normals: 3xN array of surface normals of N triangles
%   centroids__m: 3xN array of surface centroids of N triangles
%   v_rels__m_per_s: 3xN array of relative velocities of N triangles
%   deltas__rad: 1xN array of angles between the flow direction and the normals of N triangles
%   density__kg_per_m3: Scalar value of the density of the gas
%   gas_temperature__K: Scalar value of the temperature of the gas
%   surface_temperatures__K: 1xN array of the N triangles' temperatures
%   energy_accommodation_coefficients: 1xN array of the energy accommodation coefficients
%                                      of N triangles
%   particles_mass__kg: Scalar value of the mass of the particles
%   temperature_ratio_method: Scalar value of the method to calculate the temperature ratio
%                             1: Exact term according to [1]
%                             2: Hyperthermal approximation according to [1]
%                             3: Hyperthermal approximation according to [2]
%
%   Outputs:
%   aeroForce__N: 3x1 array of the aerodynamic force acting on the body in the same coordinate
%                 system as the inputs normals and centroids
%   aeroTorque__Nm: 3x1 array of the aerodynamic torque acting on the body in the same
%                   coordinate system as the inputs normals and centroids and with respect to its origin 
%
%% References:
% [1] L. H. Sentman, “Free Molecule Flow Theory and Its Application to the Determination of Aerodynamic Forces,” Defense Technical Information Center, Fort Belvoir, VA, LMSC-448514, Oct. 1961.
% [2] F. Tuttas, C. Traub, M. Pfeiffer, and W. Fichter, “Generalized Treatment of Energy Accommodation in Gas-Surface Interactions for Satellite Aerodynamics Applications,” 2024, arXiv. doi: 10.48550/ARXIV.2411.11597.
% [3] G. Koppenwallner, “Energy Accomodation Coefficient and Momentum Transfer Modeling,” HTG-TN-08-11, Dec. 2009.



%% Abbreviations
v_rels = v_rels__m_per_s;
V = vecnorm(v_rels);
rho = density__kg_per_m3;
Ti = gas_temperature__K;
Tw = surface_temperatures__K;
alpha = energy_accommodation_coefficients;
m = particles_mass__kg;

%% Constants
% Boltzmann constant
kB = 1.38064852e-23; % J K^-1

%% Local flat plate forces according to [1]
% Most probable thermal velocity of the gas
cm = sqrt(2 * kB * Ti / m);
% Molecular speed ratio 
s = V / cm;

% Intermediate values
cosdeltas = cos(deltas__rad);
scosdeltas = s .* cosdeltas;
erfcterm = erfc(-scosdeltas);
G1 = scosdeltas/sqrt(pi) .* exp(-scosdeltas.^2) + (1/2 + scosdeltas.^2) .* erfcterm;
G2 = 1/sqrt(pi) * exp(-scosdeltas.^2) + scosdeltas .* erfcterm;

% Temperature ratio
switch temperature_ratio_method
    case 1
        % Exact term according to [2]
        enum = scosdeltas .* erfcterm;
        denom = 1/sqrt(pi) * exp(-scosdeltas.^2) + enum;
        T_rat = alpha .* (2*kB*Tw)./(m*V.^2) .* s.^2 + (1-alpha) .* (1 + s.^2/2 + 1/4 * enum./denom);
    case 2
        % Hyperthermal approximation according to [2]
        T_rat = s.^2/2 .* (1 + alpha .* ((4*kB*Tw)./(m*V.^2) - 1)) + 5/4 * (1 - alpha);
    case 3
        % Hyperthermal Approximation according to [3]
        T_rat = s.^2/2 * (1 + alpha .* ((4*kB*Tw)./(m*V.^2) - 1));
    otherwise
        error('Invalid temperature ratio method');
end

% Momentum flux
p = rho/2 * cm^2 * (- ( G1 + sqrt(pi)/2 * sqrt(T_rat) .* G2 ) .* normals + s .* G2 .* (v_rels/V + cosdeltas .* normals));

%% Global forces and torques
% Forces
F = p .* areas__m2;

% Torques
tau = cross(centroids__m, F);

% Output
aeroForce__N = sum(F,2);
aeroTorque__Nm = sum(tau,2);

end