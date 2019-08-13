function [param_est , x , fc_blasso , fc_lasso , fc_lassoDual ] = FW( y , opts )

% This code is an implementation of the Frank-Wolfe algorithm [1] for the
% blasso problem.

% Inputs:
%   y : the M*1 observation vector
% opts: the parameters structure, with REQUIRED fields:
%       cplx : boolean
%       A : the M*N dictionary
%       param_grid : the 1*N vector of the discretized atom parameters subspace
%       atom : the atom generative function
%       datom : the derivative of the atom generative function
%       B : the parameters bound [ min ; max ]
%       lambda : positive real, parameter of the lasso problem
% optional fields:
%       disp : boolean for the display
%       tol :  the stopping criteria tolerence
%       maxIter : the maximum nb of iterations
%
%
%
% Ref:
%
% [1] M. Frank and P. Wolfe, 
% An algorithm for quadratic programming, Naval Research Logistics (NRL), vol.3, no. 1-2, pp. 95–110, 1956.


%% default parameters
if(~isfield(opts,'disp'))
    opts.disp=false;
end
if(~isfield(opts,'tol'))
    opts.tol=1.e-5;
end
if(~isfield(opts,'maxIter'))
    opts.maxIter=1.e2;
end

if(opts.disp)
    disp('Frank-Wolfe running...')
end

%% Initialisation
M =norm(y)^2/(2*opts.lambda);
t = M;
A = [];
param_est=[];
x=[];
residual = -y;

fc_blasso = zeros(1,opts.maxIter);
fc_lasso = zeros(1,opts.maxIter);
fc_lassoDual = zeros(1,opts.maxIter);

for iter = 1 : opts.maxIter
    
    % % % % % % % % % % % % %
    % Atom selection step % %
    % % % % % % % % % % % % %
    [ out_atom_selection ] = atom_selection( residual , opts );
    param_new = out_atom_selection.param_new;
    val_new = out_atom_selection.val;
    
    if(opts.disp)
        disp('--------')
        disp(['Iteration :',int2str(iter)])
        disp(['Selected parameter :',num2str(param_new)])
        disp(['Inner product value :',num2str(val_new)])
    end
    
    
    
    % % % % % % % % % % %
    % Convergence check %
    % % % % % % % % % % %
    if(iter>1)
        dualv = -residual/abs(val_new);
        fc_lassoDual(iter-1) = .5*norm(y)^2-.5*norm(y-opts.lambda*dualv)^2;
        dualgap = fc_lasso(iter-1)-fc_lassoDual(iter-1);
        if(dualgap<=opts.tol)
            fc_blasso = fc_blasso(1:iter-1);
            fc_lasso = fc_lasso(1:iter-1);
            fc_lassoDual = fc_lassoDual(1:iter-1);
            break;
        end
    end
    
    % % % % % % % % % % %
    % Std F.-W. update % %
    % % % % % % % % % % %
    if(opts.lambda>=abs(val_new))
        [x , t , stopflag] = std_FW_update_v1( residual , A , x , t , M , opts.lambda );
        if(stopflag)
            fc_blasso = fc_blasso(1:iter-1);
            fc_lasso = fc_lasso(1:iter-1);
            fc_lassoDual = fc_lassoDual(1:iter-1);
            break;
        end            
    else
        param_est = [ param_est , param_new ]; %#ok<AGROW>
        new_atom = opts.atom(param_new);
        [ x , t ] = std_FW_update_v2( residual , A , x , t , opts.lambda , M , val_new , new_atom , opts.cplx );
    end
    
    A = opts.atom(param_est);
    Ax  = opts.atom(param_est)*x;
    residual = Ax-y;
    fc_blasso(iter) = blasso_FObj( residual , opts.lambda , t );
    fc_lasso(iter) = lasso_FObj( residual , opts.lambda , x );
    
    if(opts.disp)
        disp(['Value of the blasso primal function: ',num2str(fc_blasso(iter))])
        disp(['Value of the lasso primal function: ',num2str(fc_lasso(iter))])
    end
    
end

if(opts.disp)
    disp('============')
    if(iter<opts.maxIter)
        disp('Convergence reached')
    else
        disp('Maximum iterations reached')
    end
    disp(['Value of the blasso primal function: ',num2str(fc_blasso(end))])
    disp(['Value of the lasso primal function: ',num2str(fc_lasso(end))])
    disp('============')
end

end


function obj = blasso_FObj( residual , lambda , t )
obj = .5*norm(residual)^2+lambda*t;
end


function obj = lasso_FObj( residual , lambda , coeff )
obj = .5*norm(residual)^2+lambda*norm(coeff,1);
end


function [ x , t , stopflag ] = std_FW_update_v1( residual , A , x , t , M , lambda )
Adiff = -A*x;
gamma = max(0,min(1,(real(-residual'*Adiff)+lambda*(t-M))/(norm(Adiff)^2)));
x = (1-gamma)*x;
t = (1-gamma)*t;
stopflag = (gamma<=1.e-10);
end

function [ x , t ] = std_FW_update_v2( residual , A , x , t , lambda ,M , Atres_new , new_atom , cplx )

if(cplx)
    coeff_new =  M*exp(1i*angle(-Atres_new));
else
    coeff_new =  M*sign(-Atres_new);
end
if(isempty(A))
    gamma = max(0,min(1,(real(-residual'*(new_atom*coeff_new))+lambda*(t-M))/(norm(new_atom*coeff_new)^2)));
    x = gamma*coeff_new;
else
    Adiff = (new_atom*coeff_new-A*x);
    gamma = max(0,min(1,(real(-residual'*Adiff)+lambda*(t-M))/(norm(Adiff)^2)));
    x = [(1-gamma)*x;gamma*coeff_new];
end
t = (1-gamma)*t+gamma*M;
end
