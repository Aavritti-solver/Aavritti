
SetFactory("OpenCASCADE");
// Define mesh size parameters
h_coarse = 0.1;
h_fine = 0.1;

// Define points with specified mesh sizes
Point(1) = {0, 0, 0, h_coarse};
Point(2) = {10.0, 0, 0, h_coarse};
Point(3) = {10.0, 1.0, 0, h_coarse};
Point(4) = {0, 1.0, 0, h_coarse};

// Define lines connecting the points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Create a curve loop and a plane surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface("domain") = {1};

// Define physical groups for boundaries and the surface
Physical Line("inlet") = {4};
Physical Line("outlet") = {2};
Physical Line("walls") = {1, 3};







