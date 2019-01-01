% ~~~~~~~~~~~~~~~~~~~~~~~~
function dYdt = twobodyEOM3D(t,Y,mu)
% ------------------------
%{
  This function calculates first and second time derivatives of r
  governed by the equation of two-body 2D motion
 
  Y    - column vector containing r and v at time t
  dYdt - column vector containing drdt and dvdt at time t

  User M-functions required: none
%}
% ~~~~~~~~~~~~~~~~~~~~~~~~

rvec = Y(1:3);
vvec = Y(4:6);

r = sqrt(rvec(1)^2+rvec(2)^2+rvec(3)^2) ;

rdotvec = vvec ;
vdotvec = -mu/r^3*rvec ;

dYdt = [rdotvec; vdotvec];
end 