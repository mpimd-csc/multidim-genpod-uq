RO = 1;                 // outer radius 
RI = 0.4;               // inner radius
AL = 2*Pi/4;            // arclength
D = 0.02*2*Pi*RI;       // mesh density
H = 0.5;                 // Extension in Z direction
DR = 0.1;             // Radius of the domain of observation

CONTDOMS = 0;           // Start of range for control domains
OBSDOMS = 10;          // Start of range for observation domains
VOLDOMS = 20;          // Start of range for physical volumes

Point(1) = {0, 0, 0, D};  // Center point
Point(2) = {0, 0, H, D};  // Center point up

// *************
// FIRST SEGMENT
// *************

Point(3) = {RI, 0, 0, D};  // East at inner radius
Point(4) = {RO, 0, 0, D};  // East at outer radius
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{3}; }}  // P5 -- inner
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{4}; }}  // P6 -- outer 

Point(7) = {RI, 0, H, D};  // East at inner radius
Point(8) = {RI+DR, 0, H, D};  // East at mid radius
Point(9) = {RO, 0, H, D};  // East at mid radius
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{7}; }}  // P10 -- inner up
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{8}; }}  // P11 -- mid up
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{9}; }}  // P12 -- outer up

// 1st SEG FLOOR

Line(1) = {3, 4};
Circle(2) = {4, 1, 6};  // outer arc
Line(3) = {6, 5};
Circle(4) = {5, 1, 3};  // inner arc

Line Loop(1) = {1, 2, 3, 4};            // The first segment
Plane Surface(1) = {1};
Physical Surface(0) = {1};              // A DUMMY FOR FENICS
Physical Surface(CONTDOMS+1) = {1};

// 1st SEG CEIL

pointpp = 6;
linepp = 4;
surfpp = 1;

// // MID RING - OBS DOMAIN
Line(1+linepp) = {1+pointpp, 2+pointpp};
Circle(2+linepp) = {2+pointpp, 2, 5+pointpp};  // middle arc
Line(3+linepp) = {5+pointpp, 4+pointpp};
Circle(4+linepp) = {4+pointpp, 2, 1+pointpp};  // inner arc

Line Loop(1+surfpp) = {1+linepp, 2+linepp, 3+linepp, 4+linepp};
Plane Surface(1+surfpp) = {1+surfpp};
Physical Surface(OBSDOMS+1) = {1+surfpp};

// // OUTER RING
Line(5+linepp) = {2+pointpp, 3+pointpp};
Circle(6+linepp) = {3+pointpp, 2, 6+pointpp};  // outer arc
Line(7+linepp) = {6+pointpp, 5+pointpp};

Line Loop(2+surfpp) = {5+linepp, 6+linepp, 7+linepp, -2-linepp};
Plane Surface(2+surfpp) = {2+surfpp};


// MID LINE FOR MESH CONTROL
pointpp = 12;
Point(pointpp+1) = {0, 0, H/2, D};  // Center point middle
Translate {0, 0, H/2} { Duplicata { Point{4}; } } // pointpp+2
Translate {0, 0, H/2} { Duplicata { Point{6}; } } // pointpp+3
Circle(8+linepp) = {2+pointpp, pointpp+1, 3+pointpp};  // outer arc middle

// 1st SEG PILLARS

linepp = linepp+8;
surfpp = surfpp+2;

Line(1+linepp) = {3, 7};
Line(2+linepp) = {4, 14};
Line(3+linepp) = {14, 9};
Line(4+linepp) = {5, 10};
Line(5+linepp) = {6, 15};
Line(6+linepp) = {15, 12};

// 1st SEG FACETS

Line Loop(1+surfpp) = {1, 14, 15, -9, -5, -13}; 
Plane Surface(1+surfpp) = {1+surfpp};  // FRONT

Line Loop(2+surfpp) = {14, 12, -17, -2}; 
Surface(2+surfpp) = {2+surfpp};  // OUTER DOWN

Line Loop(3+surfpp) = {-12, 15, 10, -18}; 
Surface(3+surfpp) = {3+surfpp};  // OUTER UP

Line Loop(4+surfpp) = {3, 16, -7, -11, -18, -17}; 
Plane Surface(4+surfpp) = {4+surfpp};  // BACK

Line Loop(5+surfpp) = {4, 13, -8, -16}; 
Surface(5+surfpp) = {5+surfpp};  // INNER

Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Volume(1) = {1};
Physical Volume(VOLDOMS+1) = {1};

// ***************
// SECOND SEGMENT
// ***************

linepp = linepp+6;
pointpp = pointpp+3;
surfpp = surfpp+5;
segid = 2;

Printf('linepp: %g', linepp);
Printf('surfpp: %g', surfpp);


// THE NEW POINTS

Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{5}; }}  // 1+pointpp -- inner bot
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{10}; }}  // 2+pointpp -- inner top 
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{11}; }}  // 3+pointpp -- mid top 
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{12}; }}  // 4+pointpp -- out top 
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{6}; }}  // 5+pointpp -- out bot 
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{pointpp}; }}  // 6+pointpp -- out mid 

// THE ARCS

Circle(1+linepp) = {5, 1, 1+pointpp};  // inner bot
Circle(2+linepp) = {10, 2, 2+pointpp};  // inner top
Circle(3+linepp) = {11, 2, 3+pointpp};  // mid top
Circle(4+linepp) = {12, 2, 4+pointpp};  // out top 
Circle(5+linepp) = {pointpp, 13, 6+pointpp};  // out mid 
Circle(6+linepp) = {6, 1, 5+pointpp};  // out bot 


// CLOSE THE FACE

Line(7+linepp) = {1+pointpp, 2+pointpp};
Line(8+linepp) = {2+pointpp, 3+pointpp};
Line(9+linepp) = {3+pointpp, 4+pointpp};
Line(10+linepp) = {4+pointpp, 6+pointpp};
Line(11+linepp) = {6+pointpp, 5+pointpp};
Line(12+linepp) = {5+pointpp, 1+pointpp};

// THE FACES
//surfpp=8
//linepp=18

Line Loop(1+surfpp) = {16, 2+linepp, -1-6-linepp, -1-linepp};   // inner
Surface(1+surfpp) = {1+surfpp};
Line Loop(2+surfpp) = {-7, 3+linepp, -2-6-linepp, -2-linepp};   // top mid ring   
Plane Surface(2+surfpp) = {2+surfpp};
Line Loop(3+surfpp) = {-11, 4+linepp, -3-6-linepp, -3-linepp};   // ceil
Plane Surface(3+surfpp) = {3+surfpp};

// linepp=18

Line Loop(4+surfpp) = {-18, 5+linepp, -4-6-linepp, -4-linepp};   // outer upper
Surface(4+surfpp) = {4+surfpp};

Line Loop(5+surfpp) = {-17, 6+linepp, -5-6-linepp, -5-linepp};   // outer lower
Surface(5+surfpp) = {5+surfpp};

Line Loop(6+surfpp) = {3, 1+linepp, -6-6-linepp, -6-linepp};   // bot
Plane Surface(6+surfpp) = {6+surfpp};

Line Loop(7+surfpp) = {7+linepp, 8+linepp, 9+linepp, 10+linepp, 11+linepp, 12+linepp}; // front
Plane Surface(7+surfpp) = {7+surfpp};

Surface Loop(segid) = {1+7, 2+7, 3+7, 4+7, 5+7, 6+7, 6};
Volume(segid) = {segid};

Physical Volume(VOLDOMS+segid) = {segid};
Physical Surface(CONTDOMS+segid) = {5+surfpp};
Physical Surface(OBSDOMS+segid) = {2+surfpp};

linepp = linepp+10;
pointpp = pointpp+5;
surfpp = surfpp+6;

//*// // **************
//*// // THIRD SEGMENT
//*// // **************
//*// 
//*// segid = 3;
//*// 
//*// Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{pointpp-4}; }}  // 1+pointpp -- inner bot
//*// Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{pointpp-3}; }}  // 2+pointpp -- inner top 
//*// Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{pointpp-2}; }}  // 3+pointpp -- mid top 
//*// Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{pointpp-1}; }}  // 4+pointpp -- out top 
//*// Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{pointpp}; }}  // 5+pointpp -- out bot 
//*// 
//*// // THE ARCS
//*// Circle(1+linepp) = {pointpp-4, 1, 1+pointpp};  // inner bot
//*// Circle(2+linepp) = {pointpp-3, 2, 2+pointpp};  // inner top
//*// Circle(3+linepp) = {pointpp-2, 2, 3+pointpp};  // mid top
//*// Circle(4+linepp) = {pointpp-1, 2, 4+pointpp};  // out top 
//*// Circle(5+linepp) = {pointpp,   1, 5+pointpp};  // out bot 
//*// 
//*// // CLOSE THE FACE
//*// 
//*// Line(6+linepp) = {1+pointpp, 2+pointpp};
//*// Line(7+linepp) = {2+pointpp, 3+pointpp};
//*// Line(8+linepp) = {3+pointpp, 4+pointpp};
//*// Line(9+linepp) = {4+pointpp, 5+pointpp};
//*// Line(10+linepp) = {5+pointpp, 1+pointpp};
//*// 
//*// // THE FACES
//*// Line Loop(1+surfpp) = {linepp-4, 2+linepp, -1-5-linepp, -1-linepp};   // inner
//*// Surface(1+surfpp) = {1+surfpp};
//*// Line Loop(2+surfpp) = {linepp-3, 3+linepp, -2-5-linepp, -2-linepp};   // top mid ring   
//*// Plane Surface(2+surfpp) = {2+surfpp};
//*// Line Loop(3+surfpp) = {linepp-2, 4+linepp, -3-5-linepp, -3-linepp};   // ceil
//*// Plane Surface(3+surfpp) = {3+surfpp};
//*// Printf('linepp: %g', linepp);
//*// Line Loop(4+surfpp) = {linepp-1, 5+linepp, -4-5-linepp, -4-linepp};   // outer
//*// Surface(4+surfpp) = {4+surfpp};
//*// Line Loop(5+surfpp) = {linepp, 1+linepp, -5-5-linepp, -5-linepp};   // bot
//*// Plane Surface(5+surfpp) = {5+surfpp};
//*// Line Loop(6+surfpp) = {6+linepp, 7+linepp, 8+linepp, 9+linepp, 10+linepp}; // front
//*// Plane Surface(6+surfpp) = {6+surfpp};
//*// 
//*// Surface Loop(segid) = {1+surfpp, 2+surfpp, 3+surfpp, 4+surfpp, 5+surfpp, 6+surfpp, surfpp};
//*// Volume(segid) = {segid};
//*// 
//*// Physical Volume(VOLDOMS+segid) = {segid};
//*// Physical Surface(CONTDOMS+segid) = {5+surfpp};
//*// Physical Surface(OBSDOMS+segid) = {2+surfpp};
//*// 
//*// linepp = linepp+10;
//*// pointpp = pointpp+5;
//*// surfpp = surfpp+6;
//*// 
//*// // **************
//*// // FINAL SEGMENT
//*// // **************
//*// 
//*// segid = 4;
//*// 
//*// // THE ARCS
//*// Circle(1+linepp) = {pointpp-4, 1, 3};  // inner bot
//*// Circle(2+linepp) = {pointpp-3, 2, 7};  // inner top
//*// Circle(3+linepp) = {pointpp-2, 2, 8};  // mid top
//*// Circle(4+linepp) = {pointpp-1, 2, 9};  // out top 
//*// Circle(5+linepp) = {pointpp,   1, 4};  // out bot 
//*// 
//*// // CLOSE THE FACE -- ITS ClOSED ALREADY
//*// 
//*// Printf('linepp: %g', linepp);
//*// Printf('surfpp: %g', surfpp);
//*// 
//*// // THE FACES
//*// Line Loop(1+surfpp) = {linepp-4, 2+linepp, -12, -1-linepp};   // inner
//*// Surface(1+surfpp) = {1+surfpp};
//*// Line Loop(2+surfpp) = {linepp-3, 3+linepp, -5, -2-linepp};   // top mid ring   
//*// Plane Surface(2+surfpp) = {2+surfpp};
//*// Line Loop(3+surfpp) = {linepp-2, 4+linepp, -9, -3-linepp};   // ceil
//*// Plane Surface(3+surfpp) = {3+surfpp};
//*// Line Loop(4+surfpp) = {linepp-1, 5+linepp, 13, -4-linepp};   // outer
//*// Surface(4+surfpp) = {4+surfpp};
//*// Line Loop(5+surfpp) = {linepp, 1+linepp, 1, -5-linepp};   // bot
//*// Plane Surface(5+surfpp) = {5+surfpp};
//*// 
//*// 
//*// Surface Loop(segid) = {1+surfpp, 2+surfpp, 3+surfpp, 4+surfpp, 5+surfpp, surfpp, 4};
//*// Volume(segid) = {segid};
//*// 
//*// Physical Volume(VOLDOMS+segid) = {segid};
//*// Physical Surface(CONTDOMS+segid) = {5+surfpp};
//*// Physical Surface(OBSDOMS+segid) = {2+surfpp};
//*// 
//*// linepp = linepp+5;
//*// 
//*// //+
//*// 
//*// Field[1] = Box;
//*// Field[1].VIn = D/2;
//*// Field[1].VOut = D;
//*// Field[1].XMax = 1;
//*// Field[1].XMin = -1;
//*// Field[1].YMax = 1;
//*// Field[1].YMin = -1;
//*// Field[1].ZMax = 0.5*D;
//*// Field[1].ZMin = -0.01;
//*// 
//*// Field[2] = Box;
//*// Field[2].VIn = D/2;
//*// Field[2].VOut = D;
//*// Field[2].XMax = RI+DR;
//*// Field[2].XMin = -RI-DR;
//*// Field[2].YMax = RI+DR;
//*// Field[2].YMin = -RI-DR;
//*// Field[2].ZMax = H;
//*// Field[2].ZMin = H-D;
//*// // 
//*// // // Field[2] = Cylinder;
//*// // // Field[2].Radius = RI+DR;
//*// // // Field[2].VIn = 0.4*D;
//*// // // Field[2].VOut = D;
//*// // // Field[2].XAxis = 0;
//*// // // Field[2].XCenter = 0;
//*// // // Field[2].YAxis = 0;
//*// // // Field[2].YCenter = 0;
//*// // // Field[2].ZAxis = 1;
//*// // // Field[2].ZCenter = 0;
//*// // //+
//*// // // Finally, let's use the minimum of all the fields as the background mesh field
//*// // Field[3] = Min;
//*// // Field[3].FieldsList = {1,2};
//*// // Background Field = 3;
//*// 
//*// Field[3] = Attractor;
//*// Field[3].NNodesByEdge = 9;
//*// Field[3].EdgesList = {10,19,29,39,2,20,30,40};
//*// 
//*// Field[4] = Threshold;
//*// Field[4].IField = 3;
//*// Field[4].LcMin = D/3;
//*// Field[4].LcMax = D;
//*// Field[4].DistMin = D/3;
//*// Field[4].DistMax = D/2;
//*// 
//*// Field[5] = Min;
//*// Field[5].FieldsList = {1,2,4};
//*// Background Field = 5;
//*// 
//*// 
//*// Point(pointpp+1) = {0, 0, H/2, D};  // Center point middle
//*// Translate {0, 0, H/2} { Duplicata { Point{4}; } } // pointpp+2
//*// Translate {0, 0, H/2} { Duplicata { Point{6}; } } // pointpp+3
//*// // Rotate { {0,0,1}, {0, 0, 0},  AL/10 } { Point{pointpp+2}; }
//*// // Rotate { {0,0,1}, {0, 0, 0},  -AL/10 } { Point{pointpp+3}; }
//*// Circle(linepp+1) = {pointpp+2, pointpp+1, pointpp+3};  // outer arc
//*// Line{linepp+1} In Volume{1};
