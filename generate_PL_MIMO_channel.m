function G = generate_PL_MIMO_channel(tx, rx, d, fc)
    lambda = physconst('LightSpeed')/(fc*1e9);
    rows = rx.Nx*rx.Ny*rx.Nz;
    cols = tx.Nx*tx.Ny*tx.Nz;
    a_tx = ones(tx.Nx*tx.Ny*tx.Nz, 1);
    a_rx = ones(rx.Nx*rx.Ny*rx.Nz, 1);
    [PL_los, ~] = path_loss(d, fc, 0, 1, 0);
    G = sqrt(PL_los)*sqrt(rows*cols)*exp(-1j*2*pi*d/lambda)*a_rx*a_tx';
end
