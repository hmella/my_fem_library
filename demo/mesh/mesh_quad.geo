// Dimensions
L  = 1.0;
H  = 1.0;
lc = 0.025;

// Exterior points
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {L, 0.0, 0.0, lc};
Point(3) = {L, H, 0.0, lc};
Point(4) = {0.0, H, 0.0, lc};

// Exterior lines
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

// Surface
Line Loop(9) = {5, 6, 7, 8};
Plane Surface(10) = {9};

Transfinite Surface{10};
Recombine Surface{10};

// Physical labels
Physical Surface(1) = {10};
Physical Line(1) = {5};
Physical Line(2) = {6};
Physical Line(3) = {7};
Physical Line(4) = {8};
