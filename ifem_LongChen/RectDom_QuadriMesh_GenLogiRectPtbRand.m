function QuadriMesh = RectDom_QuadriMesh_GenLogiRectPtbRand(...
  xa,xb,nx,yc,yd,ny,delta,status)
%% For a given rectangular domain,
% generating a quadrilateral mesh that is logically rectangular 
% and perturbation is random
% Input:  
%   xa: left  xb: right  nx: number of partitions in the x-direction 
%   yc: bottom  yd: top  ny: number of partitions in the y-direction 
%   delta: pertubation magnitude 
% Output: 
%   QuadriMesh: a structure of arrays for primary mesh info 
% JL20170111: This function is mainly for simplicity not efficiency
% James Liu, ColoState; 2012/07--2017/01

%% For status=1 (Primary mesh info) 
QuadriMesh.NumNds = (nx+1)*(ny+1);  % Number of nodes 
QuadriMesh.NumEms = nx*ny;          % Number of rectangular elements 

%%  
x = linspace(xa,xb,nx+1);
y = linspace(yc,yd,ny+1)';
[X,Y] = meshgrid(x,y);  % Labelling from bottom to top, from left to right 

%% JL20170111: TO BE REVISED FOR EFFECIENCY 
%% Perturbation in Line 30,31 
hx = (xb-xa)/nx;  hy = (yd-yc)/ny;
k = 1;
for i=1:(nx+1)
  for j=1:(ny+1)
    X(k) = X(k) + hx * sin(pi*(i-1)/nx) * delta * randn();
    Y(k) = Y(k) + hy * sin(pi*(j-1)/ny) * delta * randn();
    k = k + 1;
  end
end

%% 
QuadriMesh.xa = xa;  QuadriMesh.xb = xb;  
QuadriMesh.yc = yc;  QuadriMesh.yd = yd;
QuadriMesh.nx = nx;  QuadriMesh.hx = (xb-xa)/nx;  QuadriMesh.x = x;
QuadriMesh.ny = ny;  QuadriMesh.hy = (yd-yc)/ny;  QuadriMesh.y = y;

%% Setting the lexicongraphical order for all nodes 
QuadriMesh.node = [X(:),Y(:)];

%% Auxiliary arrays 
LblNd = reshape((1:(nx+1)*(ny+1)),ny+1,nx+1);
LblEgVert = reshape(1:(nx+1)*ny,ny,nx+1);
LblEgHori = reshape(1:nx*(ny+1),ny+1,nx);

%% Setting up element-vs-nodes (counterclockwise orientation) 
QuadriMesh.elem = zeros(QuadriMesh.NumEms,4);
lbl = LblNd(1:ny,1:nx);  % Bottom to top, left to right 
QuadriMesh.elem(:,1) = lbl(:);                         % lower-left 
QuadriMesh.elem(:,2) = QuadriMesh.elem(:,1) + (ny+1);  % lower-right 
QuadriMesh.elem(:,3) = QuadriMesh.elem(:,1) + (ny+2);  % upper-right 
QuadriMesh.elem(:,4) = QuadriMesh.elem(:,1) + 1;       % upper-left 

%% 
if (status==1)
  QuadriMesh.flag = 1;
  return;
end

%% For status=2 (Secondary mesh info)
NumEgsVert = (nx+1)*ny;
NumEgsHori = nx*(ny+1);
QuadriMesh.NumEgs = NumEgsVert + NumEgsHori;
QuadriMesh.edge = zeros(QuadriMesh.NumEgs,2);

%% Setting up vertical and horizontal edges respectively 
% Vertical edges 
lbl = LblNd(1:ny,1:nx+1);
QuadriMesh.edge(1:NumEgsVert,1) = lbl(:);
QuadriMesh.edge(1:NumEgsVert,2) = QuadriMesh.edge(1:NumEgsVert,1) + 1;
% Horizonatl edges 
lbl = LblNd(1:ny+1,1:nx);
QuadriMesh.edge(NumEgsVert+1:end,1) = lbl(:); 
QuadriMesh.edge(NumEgsVert+1:end,2) = QuadriMesh.edge(NumEgsVert+1:end,1) + (ny+1); 

%% Generating secondary mesh info on element-vs-edges 
QuadriMesh.elem2edge = zeros(QuadriMesh.NumEms,4);
lbl = LblEgVert(1:ny,1:nx);
QuadriMesh.elem2edge(:,4) = lbl(:);
QuadriMesh.elem2edge(:,2) = QuadriMesh.elem2edge(:,4) + ny;
lbl = LblEgHori(1:ny,1:nx);
QuadriMesh.elem2edge(:,1) = NumEgsVert + lbl(:);
QuadriMesh.elem2edge(:,3) = QuadriMesh.elem2edge(:,1) + 1;

%% Generating secondary mesh info on edge-vs-elements based on elem2edge 
QuadriMesh.edge2elem = zeros(QuadriMesh.NumEgs,2);
CntEmsEg = zeros(QuadriMesh.NumEgs,1);
for ie=1:QuadriMesh.NumEms
  LblEg = QuadriMesh.elem2edge(ie,1:4);
  CntEmsEg(LblEg) = CntEmsEg(LblEg) + 1;
  for k=1:4
    QuadriMesh.edge2elem(LblEg(k),CntEmsEg(LblEg(k))) = ie;
  end
end

%% Adjusting 
ig = find(QuadriMesh.edge2elem(:,1)>QuadriMesh.edge2elem(:,2));
tmp = QuadriMesh.edge2elem(ig,1);
QuadriMesh.edge2elem(ig,1) = QuadriMesh.edge2elem(ig,2);
QuadriMesh.edge2elem(ig,2) = tmp;
ig = find(QuadriMesh.edge2elem(:,1)==0);
QuadriMesh.edge2elem(ig,1) = QuadriMesh.edge2elem(ig,2);
QuadriMesh.edge2elem(ig,2) = 0;

%% Elementwise area and center 
k1 = QuadriMesh.elem(:,1);  k2 = QuadriMesh.elem(:,2);  
k3 = QuadriMesh.elem(:,3);  k4 = QuadriMesh.elem(:,4);
x1 = QuadriMesh.node(k1,1);  y1 = QuadriMesh.node(k1,2);
x2 = QuadriMesh.node(k2,1);  y2 = QuadriMesh.node(k2,2);
x3 = QuadriMesh.node(k3,1);  y3 = QuadriMesh.node(k3,2);
x4 = QuadriMesh.node(k4,1);  y4 = QuadriMesh.node(k4,2);
QuadriMesh.area = 0.5*( (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1)...
                      + (x3-x1).*(y4-y1) - (x4-x1).*(y3-y1));
QuadriMesh.EmCntr = [0.25*(x1+x2+x3+x4), 0.25*(y1+y2+y3+y4)];

%% Coefficients for the bilinear mapping 
QuadriMesh.CofA = zeros(QuadriMesh.NumEms,4);       
QuadriMesh.CofA(:,1) = x1;
QuadriMesh.CofA(:,2) = x2-x1;
QuadriMesh.CofA(:,3) = x4-x1;
QuadriMesh.CofA(:,4) = (x1+x3)-(x2+x4);
QuadriMesh.CofB = zeros(QuadriMesh.NumEms,4);
QuadriMesh.CofB(:,1) = y1;
QuadriMesh.CofB(:,2) = y2-y1;
QuadriMesh.CofB(:,3) = y4-y1;
QuadriMesh.CofB(:,4) = (y1+y3)-(y2+y4);

%% Computing elementwise diameter
diag13 = sqrt((x1-x3).^2+(y1-y3).^2);
diag24 = sqrt((x2-x4).^2+(y2-y4).^2);
QuadriMesh.diam = max(diag13,diag24);

%% Finishing secondary mesh info
QuadriMesh.flag = 2;

%% Afternotes: Mesh info is organized in the following format 
% status for the structure: 
%   1 - primary mesh info only; 
%   2 - secondary mesh info completed 
% NumNds, NumEms
% node, elem 
% NumEgs, edge, elem2edge, edge2elem
% area, EmCntr 

return;