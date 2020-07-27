function [PL_los, PL_nlos] = path_loss(d, fc, blocked, h_ut, dB)
PL_los_dB = 32.4 + 21 * log10(d) + 20 * log10(fc);
PL_nlos_dB = 22.4 + 35.3 * log10(d) + 21.3 * log10(fc) - 0.3*(h_ut-1.5);
PL_nlos_dB = max(PL_los_dB, PL_nlos_dB);
if blocked
    PL_los_dB = PL_los_dB + 20;
    PL_nlos_dB = PL_nlos_dB + 20;
end
if true(dB)
    PL_los = PL_los_dB;
    PL_nlos = PL_nlos_dB;
else
    PL_los = db2pow(-PL_los_dB);
    PL_nlos = db2pow(-PL_nlos_dB);
end
end