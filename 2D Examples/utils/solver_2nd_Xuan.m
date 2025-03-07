function solver_2nd_Xuan(Input,PF,file_name)

% Reads the model and solver parameters from the input file
[geometry,Mater,Fract,Integ] = Input();
[PHTelem,controlPts,dimBasis] = PF.Geometry(geometry);
[PHTelem,sizeBasis] = zipConforming(PHTelem,dimBasis,geometry);

scrsz = get(groot, 'ScreenSize');
hFig = figure('Position',[1 scrsz(4)/6 3*scrsz(3)/5 3*scrsz(4)/4]);
plot1 = subplot(2,2,[1,2]);
cla(plot1)
plotMesh2D(PHTelem,controlPts,geometry)
axis equal
title('Intial Mesh');

disp('Initializing boundary conditions on the initial geometry.')
[dirichlet] = PF.Boundary(PHTelem,geometry,controlPts);

disp('Precomputing shape functions and derivatives.')
[shape,dgdx,volume,gaussCord,numberElements] = cartdev(PHTelem,controlPts,geometry);

disp('History function and phase field initialization.')
Fract.constl = 2*Fract.constl;
[fenerg] = PF.History(gaussCord,Fract,numberElements,geometry); % initialize pre-existing crack
Fract.constl = Fract.constl/2;
clear gaussCord

tdisp = zeros(3*sizeBasis,1);
Integ.tfacto = 0.0;
Integ.dfacto = 0.0;
for istep = 1:Integ.nstep+1
    
    fprintf('Running... %5d/%5d\n', istep, Integ.nstep);
    if (istep < Integ.numStepsLimit)
        % Integ.tfacto = Integ.tfacto + Integ.dfacto1;
        Integ.dfacto = Integ.dfacto1;
    else
        % Integ.tfacto = Integ.tfacto + Integ.dfacto2;
        Integ.dfacto = Integ.dfacto2;
    end
    Integ.tfacto = Integ.tfacto + Integ.dfacto;
    % Begin inner iteration
    miter = 0;
    while true
        refFlag = 0;
        % initialize incremental vector tdispInc and iterative vector tdispIter
        uInc  = zeros(2*sizeBasis,1);
        [dofBC, valBC] = getDofBC(dirichlet, Integ.dfacto);
        dofFree = setdiff([1: 2*sizeBasis]', dofBC);  % Dirichlet boundary condition is only applied to displacement field
        uInc(dofBC) = valBC;

        if miter == 0
            % ==================================================
            %         SOLVING DISPLACEMENT FIELD
            % ==================================================
            % ------------ Assuming fint_{f} = 0 at beginning of each loading step
            % ------------ Using Newton-Raphson solver to solve crack phase-field
            tdispMix = [uInc(1:2*sizeBasis, 1); tdisp(2*sizeBasis+1:3*sizeBasis, 1)];
            [stiffUU, ~, elemRef] = PF.StiffUU(PHTelem,sizeBasis,numberElements,dgdx,shape,...
                    Fract,Mater,volume,tdispMix,geometry);
            KffU = stiffUU(dofFree, dofFree);
            KfdU = stiffUU(dofFree, dofBC);
            % ------------ Kff*DUf + Kfd*DUd = 0
            % ------------ DUf = [Kff]^(-1)*(-Kfd*DUd)
            uInc(dofFree) = KffU \ (-KfdU*uInc(dofBC));
            tdisp(1:2*sizeBasis) = tdisp(1:2*sizeBasis) + uInc(1:2*sizeBasis);
            % ==================================================
            %         SOLVING CRACK PHASE-FIELD
            % ==================================================
            %  ------------ Using Newton-Raphson solver to solve crack phase-field
            % fenerg = internalForces_Xuan(PHTelem,dgdx,shape,tdisp,geometry,Mater,Fract,fenerg);
            % [stiffPhiPhi,fextPhi,fintPhi] = PF.StiffPhiPhi(PHTelem,sizeBasis,numberElements,...
            %         dgdx,shape,Fract,Mater,volume,geometry,fenerg,tdisp);
            % phiInc = stiffPhiPhi\(fextPhi-fintPhi);
            % solPhi = tdisp(2*sizeBasis+1:end);
            % phiInc(solPhi > Fract.phicr) = 0.0;
            % solPhi = solPhi + phiInc;
            % solPhi(solPhi < 0.0) = 0.0;
            % solPhi(solPhi > Fract.phicr) = 1.0;
            % tdisp(2*sizeBasis+1:end) = solPhi;
            %  ------------ Using linear solver to solve crack phase-field
            fenerg = internalForces_Xuan(PHTelem,dgdx,shape,tdisp,geometry,Mater,Fract,fenerg);
            [stiffPhiPhi,fextPhi,~] = PF.StiffPhiPhi(PHTelem,sizeBasis,numberElements,...
                    dgdx,shape,Fract,Mater,volume,geometry,fenerg,tdisp);
            solPhi = stiffPhiPhi\fextPhi;
            solPhi(solPhi < 0.0) = 0.0;
            solPhi(solPhi > Fract.phicr) = 1.0;
            tdisp(2*sizeBasis+1:end) = solPhi;

            dispNorm0 = norm(uInc)/(2*sizeBasis) + 1e-10;
        else
            % ==================================================
            %         SOLVING DISPLACEMENT FIELD
            % ==================================================
            tdispMix = [uInc(1:2*sizeBasis, 1); tdisp(2*sizeBasis+1:3*sizeBasis, 1)];
            [stiffUU, fintUU, elemRef] = PF.StiffUU(PHTelem,sizeBasis,numberElements,dgdx,shape,...
                    Fract,Mater,volume,tdispMix,geometry);
            KffU = stiffUU(dofFree, dofFree);
            % ------------ Kff*dUf + Kfd*dUd + fint_U = 0 with dUd = 0
            % ------------ DUf = [Kff]^(-1)*(-fint_U)
            uIter = KffU \ (-fintUU(dofFree));
            uInc(dofFree) = uInc(dofFree) + uIter;
            tdisp(1:2*sizeBasis) = tdisp(1:2*sizeBasis) + uInc(1:2*sizeBasis);
            % ==================================================
            %         SOLVING CRACK PHASE-FIELD
            % ==================================================
            %  ------------ Using Newton-Raphson solver to solve crack phase-field
            % fenerg = internalForces_Xuan(PHTelem,dgdx,shape,tdisp,geometry,Mater,Fract,fenerg);
            % [stiffPhiPhi,fextPhi,fintPhi] = PF.StiffPhiPhi(PHTelem,sizeBasis,numberElements,...
            %         dgdx,shape,Fract,Mater,volume,geometry,fenerg,tdisp);
            % phiInc = stiffPhiPhi\(fextPhi-fintPhi);
            % solPhi = tdisp(2*sizeBasis+1:end);
            % phiInc(solPhi > Fract.phicr) = 0.0;
            % solPhi = solPhi + phiInc;
            % solPhi(solPhi < 0.0) = 0.0;
            % solPhi(solPhi > Fract.phicr) = 1.0;
            % tdisp(2*sizeBasis+1:end) = solPhi;
            %  ------------ Using linear solver to solve crack phase-field
            fenerg = internalForces_Xuan(PHTelem,dgdx,shape,tdisp,geometry,Mater,Fract,fenerg);
            [stiffPhiPhi,fextPhi,~] = PF.StiffPhiPhi(PHTelem,sizeBasis,numberElements,...
                    dgdx,shape,Fract,Mater,volume,geometry,fenerg,tdisp);
            solPhi = stiffPhiPhi\fextPhi;
            solPhi(solPhi < 0.0) = 0.0;
            solPhi(solPhi > Fract.phicr) = 1.0;
            tdisp(2*sizeBasis+1:end) = solPhi;

            errU = norm(uIter)/(2*sizeBasis) / dispNorm0;
            fprintf('    Newton-Raphson iteration %5d, err: %6.3e\n', miter, errU);
            if errU < geometry.toler
                disp('   Newton-Raphson iteration converged.');
                break;
            end
            if miter == 20
                disp('   Newton-Raphson iteration reached its maximum.');
                break;
            end
        end

        miter = miter + 1;
        clear stiffUU fintUU stiffPhiPhi fextPhi
        % [solPhiPatch] = transferFieldGlob2Loc(PHTelem,dimBasis,solPhi);
        solUPhi = reshape(tdisp, sizeBasis, 3);
        solPatch = transferFieldGlob2Loc(PHTelem,dimBasis,solUPhi);
        for iPatch = 1:geometry.numPatches
            if sum(elemRef{iPatch})>0
                % Refine and update the mesh
                refFlag = 1;
                fprintf('    In patch %5d, refining %5d elements.\n', iPatch, sum(elemRef{iPatch}))
                % [PHTelem{iPatch},controlPts{iPatch},dimBasis(iPatch),solPhiPatch{iPatch}, ...
                %     numberElements] = refineElemProjGradedIso(elemRef{iPatch},PHTelem{iPatch}, ...
                %     controlPts{iPatch},geometry,dimBasis(iPatch),solPhiPatch{iPatch},numberElements);
                [PHTelem{iPatch},controlPts{iPatch},dimBasis(iPatch),solPatch{iPatch}, ...
                    numberElements] = refineElemProjGradedIso(elemRef{iPatch},PHTelem{iPatch}, ...
                    controlPts{iPatch},geometry,dimBasis(iPatch),solPatch{iPatch},numberElements);
            end
        end

        if refFlag
            clear stiffUU
            [PHTelem,sizeBasis] = zipConforming(PHTelem,dimBasis,geometry);
            plot1 = subplot(2,2,[1,2]);
            cla(plot1)
            plotMesh2D(PHTelem,controlPts,geometry)
            axis equal
            title(['Modified Mesh for Loadstep', num2str(istep) ,' and Iteration ', num2str(miter)]);
            [dirichlet] = PF.Boundary(PHTelem,geometry,controlPts);
            [shape,dgdx,volume,gaussCord,numberElements]=cartdev(PHTelem,controlPts,geometry);
            [fenerg] = PF.History(gaussCord,Fract,numberElements,geometry); % initialize pre-existing crack
            clear gaussCord

            % solPhi = transferFieldLoc2Glob(PHTelem,sizeBasis,solPhiPatch);
            % tdisp = zeros(3*sizeBasis,1); % Solution Vector
            % tdisp(2*sizeBasis+1:end) = solPhi;
            solUPhi = transferFieldLoc2Glob(PHTelem,sizeBasis,solPatch);
            tdisp = reshape(solUPhi, 3*sizeBasis, 1);
            
        end
    end

    disp('    Print data for force-tdisp curves.')
    tic
    PF.Trac(stiffUU,tdisp(1:2*sizeBasis),Integ.tfacto,dirichlet,file_name);
    if(mod(istep,Integ.nprint) == 0)% Print results
        fprintf('Done step: %5d\n',istep);
        % plotDispPhase2D_Xuan(PHTelem,tdisp,sizeBasis,numberElements,geometry,controlPts,Mater, fenerg);
        plotDispPhase2D(PHTelem,tdisp,sizeBasis,numberElements,geometry,controlPts,Mater);
        plot1 = subplot(2,2,[1,2]);
        title(['Mesh for Loadstep ',num2str(istep),' and Iteration ',num2str(miter)]);
        saveas(hFig, ['Loadstep', num2str(istep),'.png'])
    end %if
    toc
end %istep
end