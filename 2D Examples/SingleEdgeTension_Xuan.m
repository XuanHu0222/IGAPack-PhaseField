% Script for a two dimensional plate under tension

close all
clear

addpath('./utils')
addpath('./example_data')
addpath('../nurbs/inst')

Input = @() Input_tensilePlate_Xuan;
PF.Geometry = @(geometry)initMeshTensile(geometry);
PF.Boundary = @(PHTelem,geometry,controlPts)initialBC_tensile(PHTelem,geometry);
PF.History = @(gaussCord,Fract,numberElements,geometry)history_tensile(gaussCord,...
                Fract,numberElements,geometry);
PF.Trac = @(stiffUU,tdisp,Integ,dirichlet,file_name)compTreacTension(stiffUU,tdisp,...
                Integ,dirichlet,file_name);

file_name = 'FD-tensile.txt';
output = fopen(file_name,'w');
fprintf(output,'%14.6e %14.6e\n',0,0);
fclose(output);

PF.StiffUU = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,Mater,...
    volume,tdisp,geometry)gStiffnessUU_Xuan(PHTelem,sizeBasis,...
    numberElements,dgdx,shape,Fract,Mater,volume,tdisp,geometry);
PF.StiffPhiPhi = @(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,...
    Mater,volume,geometry,fenerg,tdisp)gStiffnessPhiPhi_Xuan(PHTelem,...
    sizeBasis,numberElements,dgdx,shape,Fract,Mater,volume,geometry,fenerg);

solver_2nd_Xuan(Input,PF,file_name);