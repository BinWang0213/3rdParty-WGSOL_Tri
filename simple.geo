SetFactory("OpenCASCADE");

lc = 0.25;
Point(1) = {0,0,0, lc};
Point(2) = {3.141592653,0,0, lc};
Point(3) = {3.141592653,3.141592653,0, lc};
Point(4) = {0,3.141592653,0, lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
