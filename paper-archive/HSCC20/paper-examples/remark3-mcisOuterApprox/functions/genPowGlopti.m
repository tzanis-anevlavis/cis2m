% Generates powers of monomials UP TO degree d in the basis used by Gloptipoly 

% INPUTS:
% n - dimension
% d - degree

% OUTPUTS:
% powers - monomial powers of the Gloptipoly basis
function powers =  genPowGlopti(n,d)
    powers = [];
    for k = 0:d
        powers = [ powers ; genpow(n,k) ];
    end
end