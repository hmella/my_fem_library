// Dimensions
L  = 1.0;
H  = 1.0;
lc = 0.05;

// Exterior points
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {L, 0.0, 0.0, lc};
Point(3) = {L, H, 0.0, lc};
Point(4) = {0.0, H, 0.0, lc};

// Exterior lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surface
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Physical labels
Physical Surface(1) = {6};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
