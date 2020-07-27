function array_response = generate_array_response(N_x, N_y, N_z, delta,...
    theta, phi)
kd = 2 * pi * delta;
Nx_Ind=0:1:N_x-1;
Ny_Ind=0:1:N_y-1;
Nz_Ind=0:1:N_z-1;
Nxx_Ind=repmat(Nx_Ind,1,N_y*N_z)';
Nyy_Ind=repmat(reshape(repmat(Ny_Ind,N_x,1),1,N_x*N_y),1,N_z)';
Nzz_Ind=reshape(repmat(Nz_Ind,N_x*N_y,1),1,N_x*N_y*N_z)';
gamma_x=1j*kd*sin(theta)*cos(phi);
gamma_y=1j*kd*sin(theta)*sin(phi);
gamma_z=1j*kd*cos(theta);
gamma_comb=Nxx_Ind*gamma_x+Nyy_Ind*gamma_y + Nzz_Ind*gamma_z;
array_response=exp(gamma_comb);
end