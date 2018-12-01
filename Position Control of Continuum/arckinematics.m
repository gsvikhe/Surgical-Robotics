function [T, link] = arckinematics(c)
% ARCKINEMATICS Takes as input a vector c = [k, phi, l] and returns the corresponding
% homogenous transformation matrix T and the arc link as a
% sequence of points

k   = c(1);
phi = c(2);
l   = c(3);


T = [ cos(phi)*cos(k*l) -sin(phi) cos(phi)*sin(k*l) (cos(phi)*(1-cos(k*l)))/k;
      sin(phi)*cos(k*l)  cos(phi) sin(phi)*sin(k*l) (sin(phi)*(1-cos(k*l)))/k;
      -sin(k*l) 0 cos(k*l) sin(k*l)/k;
      0 0 0 1
    ];

link = make_link(c);
end

