hf = 0.125 * 2;
h = hf * 2;
hz = hf;
fac = 0.6;

Lf = 2;
Lxy = 5;
Lz = 4*hz;
alpha = 30;
alpha = alpha/180*Pi;

// Algorithm parameters
Mesh.Smoothing = 10; // To have a nicer mesh
Mesh.Algorithm = 5; // Delaunay for quads

// Geometry Boundary
Point(1) = {-Lxy, -Lxy, 0, h};
Point(2) = {+Lxy, -Lxy, 0, h};
Point(3) = {+Lxy, +Lf/2*Sin(alpha), 0, h};
Point(4) = {+Lxy, +Lxy, 0, h};
Point(5) = {-Lxy, +Lxy, 0, h};
Point(6) = {-Lxy, -Lf/2*Sin(alpha), 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// Crack Geometry
Point(7) = {-Lf/2*Cos(alpha), -Lf/2*Sin(alpha), 0, hf};
Point(8) = {-fac*Lf/2*Cos(alpha), -fac*Lf/2*Sin(alpha), 0, hf};
Point(9) = {+fac*Lf/2*Cos(alpha), +fac*Lf/2*Sin(alpha), 0, hf};
Point(10) = {+Lf/2*Cos(alpha), +Lf/2*Sin(alpha), 0, hf};

Point(11) = {-(1+fac)*Lf/2*Cos(alpha), -Lf/2*Sin(alpha), 0, hf};
Point(12) = {(1+fac)*Lf/2*Cos(alpha), +Lf/2*Sin(alpha), 0, hf};

Line(7) = {6, 11}; // <-- dummy divider
Line(12) = {11, 7}; // <-- dummy divider
Line(8) = {7, 8}; // <-- crack
Line(9) = {8, 9}; // <-- crack
Line(10) = {9, 10}; // <-- crack
Line(11) = {10, 12}; // <-- dummy divider
Line(13) = {12, 3}; // <-- dummy divider

//Make line loops and surfaces
Line Loop(1) = {1, 2, -13, -11, -10, -9, -8, -12, -7, 6}; // Bottom half
Line Loop(2) = {3, 4, 5, 7, 12, 8, 9, 10, 11, 13}; // Top half
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Recombine Surface{1, 2};

out[] = Extrude {0, 0, Lz} {
  // boundary
  Point{7, 10};
  // fault trace
  Line{8, 9, 10};
  Surface{1, 2};
  Layers{Lz/hz};
  Recombine;
};

// By default, the list contains the ``top'' of the extruded
// entity at index 0 and the extruded entity at index 1, followed by the
// ``sides'' of the extruded entity at indices 2, 3, etc.
// Printf("%g", out[1]);

Physical Volume("DOMAIN") = {1, 2};

Physical Surface("XMIN") = {out[16+11], out[28+4]};
Physical Surface("XMAX") = {out[16+3], out[28+2]};
Physical Surface("YMIN") = {out[16+2]};
Physical Surface("YMAX") = {out[28+3]};
Physical Surface("ZMIN") = {1, 2};
Physical Surface("ZMAX") = {out[16+0], out[28+0]};

// fault
Physical Surface("FRACTURE") = {out[4+1], out[8+1], out[12+1]}; // fault will have ID 10

Mesh.MshFileVersion = 2.2;
