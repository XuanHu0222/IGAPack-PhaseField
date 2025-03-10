function [stiffkUU, fintUU, elemRef] = gStiffnessUU_Xuan(PHTelem,sizeBasis,numberElements,dgdx,shape,Fract,Mater,volume,tdispMix,geometry)
% Note:
%   tdispMix is a mixture of tdispInc for [1:2*sizeBasis] and tdisp for [2*sizeBasis+1: 3*sizeBasis]

% Assembles the stiffness matrix and RHS (Galerkin method)

dim = geometry.dim;
nstress = geometry.nstress;
numPatches = geometry.numPatches;
ngaussX = geometry.ngaussX;
ngaussY = geometry.ngaussY;
kUU =cell(numberElements,1);
fUU =cell(numberElements,1);
elemRef = cell(numPatches,1);
elementCounter = 0;
indexCounter = 0;

for indexPatch = 1:length(PHTelem)
    
    elemRef{indexPatch} = zeros(1, length(PHTelem{indexPatch}));
    
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter =  elementCounter+1;
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1;2*sctrx],1,dim*nument);
            dispIncElmt = tdispMix(dsctrx);
            phiElmt = tdispMix(sctrx+2*sizeBasis)';
            indexCounter = indexCounter + dim^2*nument^2;
            localkUU = zeros(dim*nument,dim*nument);
            localfUU = zeros(dim*nument,1);
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    phigp = phiElmt*shape{elementCounter}(kgauss,:)';
                    % Mark element for refinement if phase field value > 0.5
                    if (phigp > geometry.threshPhi) && (PHTelem{indexPatch}(i).level < geometry.maxRefLevel)
                        elemRef{indexPatch}(i)=1;
                    end
                    
                    [Bu,~,DBu]=strainGrad(dgdx{elementCounter},nument,nstress,dim,kgauss,Mater.C);
                    
                    % Calculation of kUU
                    coeff = 1.0;
                    if phigp >= Fract.phicr
                        coeff = Fract.constk;
                    end
                    localkUU=localkUU+coeff.*(Bu'*DBu).*volume(elementCounter,kgauss);
                    localfUU=localfUU+coeff.*(Bu'*DBu*dispIncElmt).*volume(elementCounter,kgauss);
                end %jgauss
            end %igauss
            kUU{elementCounter} = localkUU;
            fUU{elementCounter} = localfUU;
        end
    end
end
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S  = zeros(1,indexCounter);
fintUU = zeros(2*sizeBasis, 1);

% Assembling the Stiffness Matrix
indexCounter = 0;
elementCounter = 0;
for indexPatch = 1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter+1;
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1;2*sctrx],1,dim*nument);            
            II(indexCounter+1:indexCounter+dim^2*nument^2) = repmat(dsctrx,1,dim*nument);
            JJ(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(repmat(dsctrx,dim*nument,1),1,dim^2*nument^2);
            S(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(kUU{elementCounter},1,dim^2*nument^2);
            indexCounter = indexCounter +dim^2*nument^2;
            fintUU(dsctrx) = fintUU(dsctrx) + fUU{elementCounter};
        end
    end
end
stiffkUU = sparse(II,JJ,S,dim*sizeBasis,dim*sizeBasis);
end