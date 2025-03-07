function [stiffkPhiPhi,fextPhi,fintPhi] = ...
    gStiffnessPhiPhi_Xuan(PHTelem,dimBasis,numberElements,dgdx,shape,Fract,Mater,volume,geometry,fenerg,tdisp)
% Assembles the stiffness matrix

dim = geometry.dim;
nstress = geometry.nstress;
kPhiPhi = cell(numberElements, 1);
fPhie = cell(numberElements, 1);
fPhii = cell(numberElements, 1); 

elementCounter = 0;
indexCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter =  elementCounter+1;
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            phiElmt = tdisp(sctrx+2*dimBasis);
            indexCounter = indexCounter + nument^2;
            localkPhiPhi = zeros(nument,nument);
            localfPhie = zeros(nument, 1);
            localfPhii = zeros(nument, 1);
            kgauss = 0;
            for ii=1:geometry.ngaussX
                for jj=1:geometry.ngaussY
                    kgauss = kgauss+1;  
                    [~,Bphi,~]=strainGrad(dgdx{elementCounter},nument,nstress,dim,kgauss,Mater.C);
                    
                    senerg = fenerg(elementCounter, kgauss);

                    % Calculation of kPhiPhi                    
                    localkPhiPhi = localkPhiPhi + Fract.cenerg*Fract.constl* (Bphi'*Bphi).*volume(elementCounter,kgauss);
                    localkPhiPhi = localkPhiPhi + (((Fract.cenerg/Fract.constl) + 0.5*senerg))...
                        .*(shape{elementCounter}(kgauss,:)'*shape{elementCounter}(kgauss,:)).*volume(elementCounter,kgauss);   
                    % calculation of fPhi
                    localfPhie = localfPhie + shape{elementCounter}(kgauss,:)'*senerg*volume(elementCounter,kgauss);
                    localfPhii = localfPhii + localkPhiPhi * phiElmt;
                end %jgaus
            end %igaus
            fPhie{elementCounter} = localfPhie;
            fPhii{elementCounter} = localfPhii;
            kPhiPhi{elementCounter} = localkPhiPhi;
        end
    end
end
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S = zeros(1, indexCounter);
fextPhi = zeros(dimBasis,1);
fintPhi = zeros(dimBasis,1);
% Assembling the Stiffness Matrix
indexCounter = 0;
elementCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter+1;
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);            
            
            II(indexCounter+1:indexCounter + nument^2) = repmat(sctrx,1,nument);
            JJ(indexCounter+1:indexCounter + nument^2) = reshape(repmat(sctrx,nument,1),1,nument^2);
            S(indexCounter+1:indexCounter + nument^2) = reshape(kPhiPhi{elementCounter},1,nument^2);
            fextPhi(sctrx) = fextPhi(sctrx) + fPhie{elementCounter};
            fintPhi(sctrx) = fintPhi(sctrx) + fPhii{elementCounter};
            indexCounter = indexCounter + nument^2;
        end
    end
end
stiffkPhiPhi = sparse(II,JJ,S,dimBasis,dimBasis);
end

