function rcs=rcs_ellipsoid(a,b,c,phi,theta)
nomi=pi*(a^2)*(b^2)*(c^2);
denomi=((a^2)*(sin(theta)^2)*(cos(phi)^2)+(b^2)*(sin(theta)^2)*(sin(phi)^2)+(c*cos(theta)^2))^2;
rcs=nomi./denomi;
end
