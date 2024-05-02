function [p_ell] = gsl_sf_legendre_Plm_e_mult(la, m, x)
% double gsl_sf_legendre_Plm_e(const int l, const int m, const double x) 
% we now accept an array of ls

if m < 0 || any(la < m) || any(x < -1.0) || any(x > 1.0)
    p_ell=0;
    error("Domain error\n");
end

nl=length(la);
if nl>1 && any(la(2:end)-la(1:end-1)<=0.0) 
    p_ell=0;
    error("Array ''la'' must be sorted in the strict accending order\n");
end
    
%/* Account for the error due to the
% * representation of 1-x.
% */
%
%/* P_m^m(x) and P_{m+1}^m(x) */

p_mm   = legendre_Pmm(m, x);
p_mmp1 = (2*m + 1) * x .* p_mm;

p_ell(nl,size(x,2))=0; % fast preallocate
%p_ell=zeros(nl,size(x,2)); % preallocate
i=1;

if la(i) == m
    p_ell(i,:) = p_mm;
    i=2;
    if i>nl
        return;
    end
end

if la(i) == m+1
    p_ell(i,:) = p_mmp1;
    i=i+1;
    if i>nl
        return;
    end
end

%/* upward recurrence: (l-m) P(l,m) = (2l-1) z P(l-1,m) - (l+m-1) P(l-2,m)
% * start at P(m,m), P(m+1,m)
% */

p_ellm2 = p_mm;
p_ellm1 = p_mmp1;

for ell=(m+2):la(end) 
    p_tmp = ((2*ell-1)*x.*p_ellm1 - (ell+m-1)*p_ellm2) ./ (ell-m);
    if la(i) == ell
        p_ell(i,:) = p_tmp;
        i=i+1;
        if i>nl
            return;
        end
    end
    p_ellm2 = p_ellm1;
    p_ellm1 = p_tmp;
end
end % function end

function p_mm = legendre_Pmm(m, x)
% double legendre_Pmm(int m, double x)
if m == 0
    p_mm = ones(1,size(x,2));
else
    p_mm = 1.0;
    %root_factor = ((1.0-x).^0.5) .* ((1.0+x).^(0.5));
    %root_factor = ((1.0-x) .* (1.0+x)).^(0.5);
    root_factor = sqrt((1.0-x) .* (1.0+x));
    fact_coeff = 1.0;
    for i=1:m
        p_mm = p_mm .* (-fact_coeff * root_factor);
        fact_coeff = fact_coeff + 2.0;
    end
end
end % function p_mm = legendre_Pmm(m, x)

