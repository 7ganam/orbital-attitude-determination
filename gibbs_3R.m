function [h e RA incl w TA a rp ra b EN]=gibbs_3R(R1,R2,R3)

deg = pi/180;
global mu

%...Data declaration for Example 5.1:
mu = 398600;
% r1 = [-294.32 4265.1 5986.7];
% r2 = [-1365.5 3637.6 6346.8];
% r3 = [-2940.3 2473.7 6555.8];
%...

%...Echo the input data to the command window:
fprintf('-----------------------------------------------------')
fprintf('\n Example 5.1: Gibbs Method\n')
fprintf('\n\n Input data:\n')
fprintf('\n  Gravitational parameter (km^3/s^2)  = %g\n', mu)
fprintf('\n  r1 (km) = [%g  %g  %g]', r1(1), r1(2), r1(3))
fprintf('\n  r2 (km) = [%g  %g  %g]', r2(1), r2(2), r2(3))
fprintf('\n  r3 (km) = [%g  %g  %g]', r3(1), r3(2), r3(3))
fprintf('\n\n');

%...Algorithm 5.1:
[v2, ierr] = gibbs(r1, r2, r3);

%...If the vectors r1, r2, r3, are not coplanar, abort:
if ierr == 1
    fprintf('\n  These vectors are not coplanar.\n\n')
    return
end

%...Algorithm 4.2:
[h e RA incl w TA a rp ra b EN] = coe_from_sv(r2,v2,mu);
% coe = [h e RA incl w TA a rp ra];


