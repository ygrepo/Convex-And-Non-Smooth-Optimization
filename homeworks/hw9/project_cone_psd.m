function B = project_cone_psd(A)
   [V,D] = eig(A);
   D(D < 0) = 0;
   B = V * D * V';
end