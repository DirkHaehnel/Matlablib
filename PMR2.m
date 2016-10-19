function [y, y1, y2] = PMR2(k1p, chi1p, chi1px, chi2p, chi2px, chip, chipx, x, flag)

if flag==1
    % alpha1p==alpha1px, eliminate k1p, alpha1p
    y1 = (chi1p + k1p * chip ^ 2 - k1p * chi1p * chi2p) ./ (1 - k1p * chi1p - k1p * chi2p - k1p.^2 * chip ^ 2 + k1p.^2 * chi1p * chi2p);
    y2 = (chi1px + x * k1p * chipx ^ 2 - x * k1p * chi1px * chi2px) ./ (1 - x * k1p * chi1px - x * k1p * chi2px - x^2 * k1p.^2 * chipx ^ 2 + x ^ 2 * k1p.^2 * chi1px * chi2px);
elseif flag==2
    % alpha2p==alpha2px, eliminate k1p, alpha2p
    y1 = (chi2p + k1p * chip ^ 2 - k1p * chi1p * chi2p) ./ (1 - k1p * chi1p - k1p * chi2p - k1p.^2 * chip ^ 2 + k1p.^2 * chi1p * chi2p);
    y2 = (chi2px + x * k1p * chipx ^ 2 - x * k1p * chi1px * chi2px) ./ (1 - x * k1p * chi1px - x * k1p * chi2px - x ^ 2 * k1p.^2 * chipx ^ 2 + x ^ 2 * k1p.^2 * chi1px * chi2px);
else
    % alphap==alphapx, eliminate k1p, alphap
    y1 = chip ./ (1 - k1p * chi1p - k1p * chi2p - k1p.^2 * chip ^ 2 + k1p.^2 * chi1p * chi2p);
    y2 = chipx ./ (1 - x * k1p * chi1px - x * k1p * chi2px - x ^ 2 * k1p.^2 * chipx ^ 2 + x ^ 2 * k1p.^2 * chi1px * chi2px);
end

y = abs(y1-y2);