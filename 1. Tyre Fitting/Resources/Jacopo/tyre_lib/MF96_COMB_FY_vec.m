% Combined Latetral force Fy
% this function remap the scalar function to its vectorial form
function [fy_vec] = MF96_COMB_FY_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

    fy_vec = zeros(size(alpha_vec));
    for i = 1:length(alpha_vec)
        % Precode
        fy0 = MF96_FY0(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
        [~, Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
        % Main Code
        fy_vec(i) = fy0 * Gyk + SVyk;
    end

end
