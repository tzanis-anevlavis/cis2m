% Returns a list of monomials up to degree d in sdpvariable x indexed in the basis used by Gloptipoly

% INPUTS: x - sdpvariable
%         d - maximum total degree of the monomials

% OUTPUTS: mlist - list of monomials
function mlist = monolistYalToGlop(x,d)
    n = numel(x);
    mlist1 = repmat(x',numMonTot(n,d),1).^genPowGlopti(n,d);
    mlist = mlist1(:,1);
    for i = 2:n
        mlist = mlist.*mlist1(:,i);
    end
end