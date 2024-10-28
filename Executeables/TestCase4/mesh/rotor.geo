// Enable OpenCASCADE geometry kernel
SetFactory("OpenCASCADE");
Mesh.ScalingFactor = 1e-3;
//+
meshFactor=0.6;
// Merging the rotor geometry
Merge "rotor.brep";

// Define curve loops for the two surfaces
Curve Loop(1) = {3, -2};
Plane Surface(1) = {1}; // Surface 1

Curve Loop(2) = {1};
Plane Surface(2) = {2}; // Surface 2 (to be subtracted from Surface 1)

// Subtract Surface 2 from Surface 1
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }

//+
Physical Curve("air") = {2};
//+
Physical Curve("source") = {1};
//+
Physical Curve("wall") = {3};
//+
Physical Surface("Domain") = {1};

//+
Transfinite Curve {1} = 100*meshFactor Using Progression 1;
Transfinite Curve {2,3} = 100*meshFactor Using Progression 1;
