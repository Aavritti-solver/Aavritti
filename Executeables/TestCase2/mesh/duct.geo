Merge "duct.brep";
//+
Mesh.ScalingFactor = 1e-3;
//+
meshFactor=0.6;
//+
SetFactory("OpenCASCADE");
//+
Physical Surface("inlet") = {5};
//+
Physical Surface("outlet") = {13};
//+
Physical Surface("walls") = {1,2,3,4,6,7,8,9,10,11,12};
//+
Physical Volume("Domain") = {1};
// Meshing
//+
Transfinite Curve {12,15} = 40*meshFactor Using Progression 1;
Transfinite Curve {14} = 30*meshFactor Using Progression 1;
//+
Transfinite Curve {3,6,9,11} = 24*meshFactor Using Progression 1;
Transfinite Curve {18,19,17,16} = 24*meshFactor Using Progression 1;
//+
Transfinite Curve {5,8,2,1} = 48*meshFactor Using Progression 1;
Transfinite Curve {25,23,20,21} = 48*meshFactor Using Progression 1;
//+
Transfinite Curve {4,10,13,7} = 24*meshFactor Using Progression 1;
Transfinite Curve {24,27,22,26} = 24*meshFactor Using Progression 1;
