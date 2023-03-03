// Gmsh

depth = 0.150;
length = 0.008;
height = 0.070;

x0 = depth/2;
y0 = length/2;
z0 = height/2;

ofs = 0;

//scaling parameters
cc1=1;
dc=cc1*30*100;
cc=8;
sc=1;

//grid size
nx =    x0*dc*sc;
ny = cc*y0*dc*sc;
nz =    z0*dc*sc;

bC=1/2;
bumpCx=cc1*bC/nx*sc^2;
bumpCy=cc*bC/ny*sc^2;
bumpCz=cc1*bC/nz*sc^2;

Nx=nx+1;
Ny=ny+1;
Nz=nz+1;
cl3 = 2*x0/Nx;

//geometry
Point(1) = { -x0+ofs, -y0, -z0, cl3};
Point(2) = { -x0+ofs, y0, -z0, cl3};
Point(3) = { -x0+ofs, y0,  z0, cl3};
Point(4) = { -x0+ofs, -y0,  z0, cl3};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

ext1[] = Extrude {2*x0, 0, 0} { Surface{1};};
//x--Depth
Transfinite Line {11,12,16,20} = Nx Using Bump bumpCx;
//y--Vertical
Transfinite Line {1,3,6,8} = Ny Using Bump bumpCy;
//z--Horizontal
Transfinite Line {2,4,7,9} = Nz Using Bump bumpCz;

Transfinite Volume "*";
Recombine Volume "*";
Transfinite Surface "*";
Recombine Surface "*";

Physical Volume("fullbox") = {ext1[1]};
//boundaries
Physical Surface("front") = {1};
Physical Surface("back") = {ext1[0]};
Physical Surface("left") = {ext1[2]};
Physical Surface("top") = {ext1[3]};
Physical Surface("right") = {ext1[4]};
Physical Surface("bot") = {ext1[5]};
