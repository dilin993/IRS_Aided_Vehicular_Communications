function G = generate_MIMO_channel(tx, rx, d,fc,blocked,beta, h_ut, PL)
    lambda = physconst('LightSpeed')/(fc*1e9);
    rows = rx.Nx*rx.Ny*rx.Nz;
    cols = tx.Nx*tx.Ny*tx.Nz;
    a_tx = generate_array_response(tx.Nx, tx.Ny, tx.Nz, tx.delta,...
        tx.departure_theta,tx.departure_phi);
    a_rx = generate_array_response(rx.Nx, rx.Ny, rx.Nz, rx.delta,...
        rx.arrival_theta,rx.arrival_phi);
    [PL_los, PL_nlos] = path_loss(d, fc, blocked, h_ut, 0);
    Glos = exp(-1j*2*pi*d/lambda)*a_rx*a_tx';
    Gnlos =(1/sqrt(2))*(randn(rows,cols) + ...
        1i*randn(rows,cols));
    if true(PL)
        Glos = sqrt(PL_los)*Glos;
        Gnlos =  sqrt(PL_nlos)*Gnlos;
    end
    G = (1/sqrt(beta+1))*Glos + sqrt(beta/(1+beta))*Gnlos;
    G = sqrt(db2pow(rx.gain)*db2pow(tx.gain))*G;
end





