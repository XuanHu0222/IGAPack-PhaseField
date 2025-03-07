function fenerg = internalForces_Xuan(PHTelem,dgdx,shape,tdisp,geometry,Mater,Fract,fenerg)
% Compute the circumferential stress

dim = geometry.dim;
nstress = geometry.nstress;

elementCounter = 0;
for indexPatch = 1:length(PHTelem)
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            elementCounter =  elementCounter+1;
            nument = size(PHTelem{indexPatch}(i).C,1);
            sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1;2*sctrx],1,dim*nument);
            dispElmt = tdisp(dsctrx);
            sizeBasis = numel(tdisp) / 3;
            phiElmt = tdisp(2*sizeBasis+sctrx);

            % Calculate strains and stresses at integration points
            kgauss=0;
            for ii=1:geometry.ngaussX
                for jj=1:geometry.ngaussY
                    kgauss=kgauss+1;
                    
                    phigp = phiElmt'*shape{elementCounter}(kgauss,:)';
                    dforce = 0.0;
                    if (phigp >= exp(-Fract.alpha)) && (phigp <= Fract.phicr) 
                        [Bu,Bphi,~]=strainGrad(dgdx{elementCounter},nument,nstress,dim,kgauss,Mater.C);
                        strainElmt = Bu*dispElmt;
                        stressElmt = Mater.C*strainElmt;
                        gradphiElmt = Bphi*phiElmt;
                        n1 = gradphiElmt(1) / (norm(gradphiElmt) + 1.0e-10);
                        n2 = gradphiElmt(2) / (norm(gradphiElmt) + 1.0e-10);
                        stre_tt = stressElmt(1)*n2^2 - 2.0*stressElmt(3)*n1*n2 + stressElmt(2)*n1^2;
                        stre_tt = 0.5 * (stre_tt + abs(stre_tt));
                        dforce = 2.0*pi / Mater.E_frac * stre_tt^2;
                    end
                    dforce = max(dforce, fenerg(elementCounter,kgauss));
                    fenerg(elementCounter,kgauss) = dforce;
                end %igaus
            end %jgaus
        end
    end
end

end %endfunction