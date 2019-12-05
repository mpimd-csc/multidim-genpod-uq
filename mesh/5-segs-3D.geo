RO = 1;                 // outer radius 
RI = 0.4;               // inner radius
AL = 2*Pi/4;            // arclength
D = 0.01*2*Pi*RI;       // mesh density
H = 1;                 // Extension in Z direction
ROBSD = 0.1;             // Radius of the domain of observation

Point(1) = {0, 0, 0, D};  // Center point

// FIRST SEGMENT

Point(2) = {RI, 0, 0, D};  // East at inner radius
Point(3) = {RO, 0, 0, D};  // East at outer radius
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{2}; }}  // P4 -- inner
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{3}; }}  // P5 -- outer 

Line(1) = {2, 3};
Circle(2) = {3, 1, 5};  // outer arc
Line(3) = {5, 4};
Circle(4) = {4, 1, 2};  // inner arc

Line Loop(1) = {1, 2, 3, 4}; // The first segment
Plane Surface(1) = {1};
Physical Surface(0) = {1};

// SECOND SEGMENT
lastline = 4;
pointpp = 5;
segid=1;
surfpp = 1;

Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{pointpp-1}; }}  // P4+pointpp -- inner
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{pointpp}; }}  // P5+pointpp -- outer 

Circle(lastline+1) = {pointpp, 1, 2+pointpp};  // outer arc
Line(lastline+2) = {2+pointpp, 1+pointpp};
Circle(lastline+3) = {1+pointpp, 1, pointpp-1};  // inner arc

Line Loop(1+surfpp) = {-lastline+1, lastline+1, lastline+2, lastline+3}; // The first segment
Plane Surface(1+surfpp) = {1+surfpp};
Physical Surface(0+surfpp) = {1+surfpp};

// THIRD SEGMENT
lastline = lastline+3;
pointpp = pointpp+2;
surfpp = surfpp+1;

Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{pointpp-1}; }}  // P4+pointpp -- inner
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{pointpp}; }}  // P5+pointpp -- outer 

Circle(lastline+1) = {pointpp, 1, 2+pointpp};  // outer arc
Line(lastline+2) = {2+pointpp, 1+pointpp};
Circle(lastline+3) = {1+pointpp, 1, pointpp-1};  // inner arc

Line Loop(1+surfpp) = {-lastline+1, lastline+1, lastline+2, lastline+3}; // The first segment
Plane Surface(1+surfpp) = {1+surfpp};
Physical Surface(0+surfpp) = {1+surfpp};

// // LAST SEGMENT

lastline = lastline+3;
surfpp = surfpp+1;

Circle(lastline+1) = {2+pointpp, 1, 3};  // outer arc
Circle(lastline+2) = {2, 1, 1+pointpp};  // inner arc

Line Loop(1+surfpp) = {-lastline+1, lastline+1, -1, lastline+2}; // The first segment
Plane Surface(1+surfpp) = {1+surfpp};
Physical Surface(0+surfpp) = {1+surfpp};


// // THE INNER REGION

Line Loop(2+surfpp) = {4, 7, 10, 12};
Plane Surface(2+surfpp) = {2+surfpp};
Physical Surface(1+surfpp) = {2+surfpp};


// // THE OBSERVATION DOMAIN

lastline = lastline+2;
surfpp = surfpp+2;
pointpp = pointpp+2;

Point(pointpp+1) = {0, 0, H, D};  // Center point at top
Point(pointpp+2) = {ROBSD, 0, H, D};
Point(pointpp+3) = {0, ROBSD, H, D};
Point(pointpp+4) = {-ROBSD, 0, H, D};
Point(pointpp+5) = {0, -ROBSD, H, D};

Circle(lastline+1) = {pointpp+2, pointpp+1, pointpp+3};
Circle(lastline+2) = {pointpp+3, pointpp+1, pointpp+4};
Circle(lastline+3) = {pointpp+4, pointpp+1, pointpp+5};
Circle(lastline+4) = {pointpp+5, pointpp+1, pointpp+2};

Line Loop(1+surfpp) = {lastline+1, lastline+2, lastline+3, lastline+4};
Plane Surface(1+surfpp) = {1+surfpp};
Physical Surface(0+surfpp) = {1+surfpp};


// // // THE OUTER BOUNDARY
// Physical Line(5) = {2, 5, 8, 11};


// GOING 3D

Extrude {0, 0, H} {
  Surface{1}; 
  Surface{2}; 
  Surface{3}; 
  Surface{4}; 
  Surface{5}; 
}
Physical Volume(30) = {1};
Physical Volume(31) = {2};
Physical Volume(32) = {3};
Physical Volume(33) = {4};
Physical Volume(34) = {5};
//+
// SetFactory("OpenCASCADE");
// Circle(84) = {0, 0, H, ROBSD, 0, 2*Pi};
// //+
// Line Loop(6) = {84};
// //+
// Line Loop(7) = {84};
// //+
// Plane Surface(123) = {6, 7};
// Physical Surface(21) = {123};
