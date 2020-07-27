function d = dist_3d(p1, p2)
d = sqrt(sum((p1 - p2).^2, 1));
end