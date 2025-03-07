function [dofBC, valBC] = getDofBC(dirichlet,dfacto)
% get the dofs and their corresponding values on dirichlet boundary
dim = 2;
fixedPts = length(dirichlet.XY);
dofBC = [];
valBC = [];

for ivfix = 1:fixedPts
    lnode = dirichlet.XY(ivfix);    
    for idofn =1:dim
        if(dirichlet.ValXY(ivfix,idofn) == 1)
            itotv =(lnode-1)*dim +idofn;
            dispFull = dirichlet.restrainedPts(ivfix,idofn).*dfacto;
            dofBC = [dofBC, itotv];
            valBC = [valBC, dispFull];
        end
    end
end
end %endfunction