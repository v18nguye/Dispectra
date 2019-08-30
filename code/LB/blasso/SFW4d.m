function [param_est , x , fc_blasso , fc_lasso , fc_lassoDual ] = SFW4d( y , opts)

% This code is an implementation of the Sliding-Frank-Wolfe algorithm [1] for the
% 4 dimensional blasso problem i.e. these are four estimated parameters.

% Inputs:
%   y : the M*1 observation vector
% opts: the parameters structure, with REQUIRED fields:
%       cplx : boolean
%       A : the M*N dictionary
%       param_grid : the 4*N vector of the discretized atom parameters subspace
%       atom : the atom generative function
%       datom : the derivative of the atom generative function
%       B : the parameters bound [ min ; max ]
%       lambda : positive real, parameter of the lasso problem
%       mergeStep : real, vicinity of candidate atoms to be merged
% optional fields:
%       disp : boolean for the display
%       tol :  the stopping criteria tolerence
%       maxIter : the maximum nb of iterations
%
%
%
% Ref:
%
% [1] Q. Denoyelle, V. Duval, G. Peyré, E. Soubies,
% The sliding frank-wolfe algorithm and its application to super-resolution microscopy. arXiv preprint arXiv:1811.06416.

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
    disp('Sliding-Frank-Wolfe running...')
end

addpath(genpath('./utils'))

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

opts_fista.lambda = opts.lambda;
opts_fista.maxIter = 1000;
opts_fista.tol = 1.e-6;
opts_fista.disp = false;


for iter = 1 : opts.maxIter
    
    % % % % % % % % % % % % %
    % Atom selection step % %
    % % % % % % % % % % % % %
    % add a supplementary theta.
    [ out_atom_selection ] = atom_selection_4d( residual , opts );
    param_new = out_atom_selection.param_new;
    val_new = out_atom_selection.val;
    
    if(opts.disp)
        disp('--------')
        disp(['Iteration :',int2str(iter)])
        disp(['x1: ',num2str(param_new(1,1)),' x2: ',num2str(param_new(2,1))])
        disp(['x3: ',num2str(param_new(3,1)),' x4: ',num2str(param_new(4,1))])
        disp(['Inner product value :',num2str(val_new)])
    end
    
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
    % modify alphas
    if(opts.lambda>=abs(val_new))
        [ x , ~ , stopflag] = std_FW_update_v1( residual , A , x , t , M , opts.lambda );
        if(stopflag)
            fc_blasso = fc_blasso(1:iter-1);
            fc_lasso = fc_lasso(1:iter-1);
            fc_lassoDual = fc_lassoDual(1:iter-1);
            break;
        end            
    else
        param_est = [ param_est , param_new ]; %#ok<AGROW>
        new_atom = opts.atom(param_new);
        [ x , ~ ] = std_FW_update_v2( residual , A , x , t , opts.lambda , M , val_new , new_atom , opts.cplx );
    end
    A = opts.atom(param_est);
    
    
    % % % % % % % % % %
    % Fista  update % %
    % % % % % % % % % %
    %update alphas
    opts_fista.L = max(eig(A'*A));
    opts_fista.xinit = x;
    opts_fista.A=A;
    x = fista(y,opts_fista);
    
    
    % % % % % % % % % %
    % Joint  update % %
    % % % % % % % % % %
    % update alpha and theta together
    [ param_est , x , ~ ] = Joint_updt( y , param_est , x , opts.lambda , opts.B , opts.atom , opts.datom , opts.cplx );
    
    
    % % % % %
    % Merge %
    % % % % %
    % merge many spikes that are closed 
    [ param_est , x , t ] = merge( y , param_est , x , opts.lambda , M , opts.atom , opts.datom , opts.B , opts.mergeStep , opts.cplx );
    
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


function [ param , coeff , t ] = Joint_updt( y , param , coeff , lambda , B , atom , datom , cplx )

% the function's arguments:
%
% - param = shape(4,k)
%       + 4: four parameters x1, x2, x3, x4. Each row corresponds to a
%       range of the parameter's values.
%       + k: a number of combinations of the parameters arranged in a column.
%
%               " A*coeff = y "
% - coeff = shape(k,1)
%       + k: a number of combinations of the parameters arranged in a
%       column.
%
% - lamda = shape(1,1)
%       + min || y - [atom(theta_1)*coeff_1 + ... atom(theta_i)*coeff_i + ... atom(theta_k)*coeff_k] ||^2 + lamda*coeff_i + ...
%            where i = 1 ... k

n = length(param(1,:)); % the number of parameters's combinations.

% - vi = shape(3*k,1) = [vi_1 vi_2 .. vi_k abs(coeff_vi_1) abs(coeff_vi_2) ..
% abs(coeff_vi_k) angle(coeff_vi_1) angle(coeff_vi_2) ..
% angle(coeff_vk_k)].T
%
% where:
%   + k is the number of parameters's combinations
%   + i is the ith variables i = {1,2,3,4}

v1 = [param(1,:)';abs(coeff);angle(coeff)]; % for the first parameter.
v2 = [param(2,:)';abs(coeff);angle(coeff)]; % ...
v3 = [param(3,:)';abs(coeff);angle(coeff)];
v4 = [param(4,:)';abs(coeff);angle(coeff)];

% v = shape(3*k,4)
%
v = [v1 v2 v3 v4]; % Concatenation each column vi


A1 = [ -eye(n) , zeros(n,2*n), zeros(n,9*n) ;... % for the first parameter.
    eye(n) , zeros(n,2*n), zeros(n,9*n) ;...
    zeros(n) , -eye(n) , zeros(n,n), zeros(n, 9*n)];
b1 = [ -ones(n,1)*B(1,1,1) ; ones(n,1)*B(2,1,1) ; zeros(n,1) ];

A2 = [ zeros(n,3*n),-eye(n) , zeros(n,2*n),zeros(n,6*n)  ;... % ...
    zeros(n,3*n), eye(n) , zeros(n,2*n), zeros(n,6*n);...
    zeros(n,3*n), zeros(n) , -eye(n) , zeros(n,n), zeros(n,6*n)];
b2 = [ -ones(n,1)*B(1,2,1) ; ones(n,1)*B(2,2,1) ; zeros(n,1) ];

A3 = [ zeros(n,6*n),-eye(n) , zeros(n,2*n),zeros(n,3*n)  ;...
    zeros(n,6*n), eye(n) , zeros(n,2*n), zeros(n,3*n);...
    zeros(n,6*n), zeros(n) , -eye(n) , zeros(n,n), zeros(n,3*n)];
b3 = [ -ones(n,1)*B(1,1,2) ; ones(n,1)*B(2,1,2) ; zeros(n,1) ];

A4 = [ zeros(n,9*n) ,-eye(n) , zeros(n,2*n);...
     zeros(n,9*n),eye(n) , zeros(n,2*n);...
    zeros(n,9*n) , zeros(n) , -eye(n) , zeros(n,n)];
b4 = [ -ones(n,1)*B(1,2,2) ; ones(n,1)*B(2,2,2) ; zeros(n,1) ];


A = [A1; A2; A3; A4];

b = [b1; b2; b3; b4];

fobj = @(v0) joint_cost( y , v0 , lambda , atom , datom , cplx );
options = optimoptions(@fmincon,'Display','off','GradObj','on','DerivativeCheck','off','Algorithm','sqp');
[ v ] = fmincon(fobj,v,A,b,[],[],[],[],[],options);

param = v(1:n,:)';
if(cplx)
    coeff = v(n+1:2*n,:).*exp(1i*v(2*n+1:3*n,:));
else
    coeff = real(v(n+1:2*n,1).*exp(1i*v(2*n+1:3*n,1)));
end
t = norm(coeff,1);

end


function [fc,grad] = joint_cost( y , v , lambda , atom , datom , cplx )

l = size(v);
l = l(1)/3;

theta = (v(1:l,:))';
alpha = v(l+1:2*l,1);
gamma = v(2*l+1:3*l,1);

coeff = alpha.*exp(1i*gamma);

A = atom(theta);
dA = datom(theta);
res = y-A*coeff;
fc = .5*norm(res,2)^2+lambda*norm(alpha,1);

coeff_4d = [];

for i = 1: length(coeff)
   
    coeff_4d = [coeff_4d coeff(i,1)*ones(1,4)];
    
end

grad_theta = -real(res'*dA*diag(coeff_4d));
grad_alpha = -(real(res'*A*diag(exp(1i*gamma)))).'+lambda;

if(cplx)
    grad_gamma = imag(res'*A*diag(coeff)).';
else
    grad_gamma = zeros(l,1);
end

grad = [grad_theta(1:4:length(grad_theta))';grad_alpha;grad_gamma;...
    grad_theta(2:4:length(grad_theta))';grad_alpha;grad_gamma;...
    grad_theta(3:4:length(grad_theta))';grad_alpha;grad_gamma;...
    grad_theta(4:4:length(grad_theta))';grad_alpha;grad_gamma;...
    ];

end


function [ param , coeff , t ] = merge( y , param , coeff , lambda , M , atom , datom , B , dist_step , cplx )
t = norm(coeff,1);

opts_fista.lambda = lambda;
opts_fista.maxIter = 1000;
opts_fista.tol = 1.e-8;
opts_fista.disp = false;

while true
    
    residual = atom(param)*coeff-y;
    InitObj = lasso_FObj( residual , lambda , coeff );
    Ddist = triu(squareform(pdist(param')));
    Ddist(Ddist==0)=nan;
    Ddist(Ddist>dist_step) = nan;
    
    
    while (nnz(~isnan(Ddist(:,:)))~=0)
        
        m = min(Ddist(:));
        [p1,p2] = find(Ddist==m,1);
        idx_p = sort([p1,p2]);
        
        param_new = mean(param(:,idx_p),2);
                
        param_temp = param;
        param_temp(:,idx_p(2))=[];
        param_temp(:,idx_p(1))=[];
        coeff_temp = coeff;
        coeff_temp(idx_p(2))=[];
        coeff_temp(idx_p(1))=[];
        t_temp = norm(coeff_temp,1);
        
              
        A_temp = atom(param_temp);
        residual = A_temp*coeff_temp-y;
        
        % atom selection
        fObj = @(param) min_scal_prod(residual,param,atom,datom);
        A = [ -1  0 0 0;1 0 0 0 ; 0 -1 0 0;0 1 0 0;0 0 -1 0; 0 0 1 0; 0 0 0 -1; 0 0 0 1];
        b = [-B(1,1,1); B(2,1,1) ;-B(1,2,1); B(2,2,1); -B(1,1,2); B(2,1,2) ;-B(1,2,2); B(2,2,2)];
        options = optimoptions(@fmincon,'Display','off','GradObj','on','Algorithm','sqp');
        param_new = fmincon(fObj,param_new,A,b,[],[],[],[],[],options);
        val_new = scal_prod( residual , param_new , atom );
        if(lambda>=abs(val_new))
            [ coeff_temp , ~] = std_FW_update_v1( residual , A_temp , coeff_temp , t , M , lambda );
        else
            param_temp = [ param_temp , param_new ]; %#ok<AGROW>
            new_atom = atom(param_new);
            [ coeff_temp , ~ ] = std_FW_update_v2( residual , A_temp , coeff_temp , t , lambda , M , val_new , new_atom , cplx );
        end
        A_temp = atom(param_temp);
        
        % coeff update (LASSO)
        opts_fista.L = max(eig(A_temp'*A_temp));
        opts_fista.xinit = coeff_temp;
        opts_fista.A=A_temp;
        coeff_temp = fista(y,opts_fista);
        
        % test
        Fc_new = lasso_FObj( A_temp*coeff_temp-y , lambda , coeff );
        
        if( InitObj - Fc_new >=0)
            disp('Merge')
            param = param_temp;
            coeff = coeff_temp;
            t = t_temp;
            break
        else
            Ddist(idx_p(1),idx_p(2)) = nan;
        end
    end
    if(nnz(~isnan(Ddist))==0)
        break
    end
end

end
