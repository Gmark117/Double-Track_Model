% Combined Longitudinal force Fx
% this function remap the scalar function to its vectorial form
function [fx_vec] = MF96_COMB_FX_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

    fx_vec = zeros(size(kappa_vec));
    for i = 1:length(kappa_vec)
        % precode
        fx0 = MF96_FX0(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
        [Gxa, ~, ~] = MF96_FXFYCOMB_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
        % main code
        fx_vec(i) = fx0 * Gxa;
    end

end