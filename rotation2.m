%function to determine the elements of rotation matrix
%system: omega phi kappa
function [r] = rotation2 (omega1,phi1,kappa1)

omega1=deg2rad(omega1);
phi1=deg2rad(phi1);
kappa1=deg2rad(kappa1);

r=[cos(phi1)*cos(kappa1) -cos(phi1)*sin(kappa1) sin(phi1);
    cos(omega1)*sin(kappa1)+sin(omega1)*sin(phi1)*cos(kappa1) cos(omega1)*cos(kappa1)-sin(omega1)*sin(phi1)*sin(kappa1) -sin(omega1)*cos(phi1);
    sin(omega1)*sin(kappa1)-cos(omega1)*sin(phi1)*cos(kappa1) sin(omega1)*cos(kappa1)+cos(omega1)*sin(phi1)*sin(kappa1) cos(omega1)*cos(phi1)];

r=transpose(r);