function fe = electronScatteringFactor(Z,q)
% ELECTRONSCATTERINGFACTOR Retrieve the parameterized electron scattering
% factor
% Z = atomic number
% q = scattering vector magnitude (1/d)
% 
% Parameterization from Kirkland is arrays given as
% a1 b1 a2 b2
% a3 b3 c1 d1
% c2 d2 c3 d3;
% fe(q) = sum_i_N (ai/(q^2 + bi) + sum_i_N ci*exp(-di*q^2);
% Parameters are in units of Angstroms

load('fparams.mat','fparams')
fparamsEl = fparams(Z,:)'; 
a = [1 3 5];
b = [2 4 6];
c = [7 9 11];
d = [8 10 12];
fe = zeros(numel(q),1);
for ii = 1:3
    fe(:) = fe ...
        + fparamsEl(a(ii))./(q(:).^2 + fparamsEl(b(ii))) ...
        + fparamsEl(c(ii)).*exp(-fparamsEl(d(ii)).*q(:).^2);
end

end