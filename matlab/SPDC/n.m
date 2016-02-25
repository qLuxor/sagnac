function out = n( l,axis,T )
%N Summary of this function goes here
%   Detailed explanation goes here

    function out = ktp( l,axis )
        l = l*1e6; % transform wavelength to micron
        if axis == 'x'
            out = sqrt(3.291 + 0.0414./(l.^2 - 0.03978) + 9.35522./(l.^2 - 31.45571));
        elseif axis == 'y'
            out = sqrt(3.45018 + 0.04341./(l.^2 - 0.04597) + 16.98825./(l.^2 - 39.43799));
        elseif axis == 'z'
            out = sqrt(4.59423 + 0.06206./(l.^2 - 0.04763) + 110.80672./(l.^2 - 86.12171));
        end
    end

    function out = T_corr_ktp(wl,axis,T)
        if axis=='y'
            n_1 = 6.2897e-6 + 6.3061e-6./(wl*1e6) - 6.0629e-6./((wl.*1e6).^2) + 2.6486e-6./((wl.*1e6).^3);
            n_2 = -0.14445e-8 + 2.2244e-8./(wl.*1e6) - 3.5770e-8./((wl.*1e6).^2) + 1.3470./((wl.*1e6).^3);
            out = n_1 .* (T-20) + n_2 .* (T-20).^2;
        elseif axis=='z'
            n_1 = 9.9587e-6 + 9.9228e-6./(wl.*1e6) - 8.9603e-6./((wl.*1e6).^2) + 4.1010e-6./((wl.*1e6).^3);
            n_2 = -1.1882e-8 + 10.459e-8./(wl.*1e6) - 9.8138e-8./((wl.*1e6).^2) + 3.1481./((wl.*1e6).^3);
            out = n_1 .* (T-20) + n_2 .* (T-20).^2;
        else
            out = 0;
        end
    end

    if nargin < 3
        T = 20;
    end

    out = ktp(l,axis) + T_corr_ktp(l,axis,T);

end

