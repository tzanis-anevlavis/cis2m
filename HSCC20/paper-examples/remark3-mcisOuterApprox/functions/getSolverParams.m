% Returns solver parameters as a Yalmip sdpsettings structure
% Parameter selected only for MOSEK and SeDuMi, otherwise using Yalmip's
% default ones

function options = getSolverParams( SDPsolver )

if(strcmpi(SDPsolver,'mosek'))
    if(~exist('mosekopt.m'))
        disp('***MOSEK not found, leaving the solver selection to Yalmip***')
        SDPsolver = '';
    end
end

if(strcmpi(SDPsolver,'sedumi'))
    if(~exist('sedumi.m'))
        disp('***SeDuMi not found, leaving the solver selection to Yalmip***')
        SDPsolver = '';
    end
end

switch SDPsolver
    case 'mosek'
        disp('SDP Solver: MOSEK')
        options = sdpsettings('solver','mosek-sdp',...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-12,...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS',1e-12,...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-12);
    case 'sedumi'
        disp('SDP Solver: SeDuMi')
        options = sdpsettings('solver','sedumi','sedumi.eps',1e-12);
    case ''
        options = [];
    otherwise
        disp(['SDP Solver: ' SDPsolver])
        options = sdpsettings('solver',SDPsolver);
end


