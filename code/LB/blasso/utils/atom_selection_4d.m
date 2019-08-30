function [ out ] = atom_selection_4d( residual , opts )

% The code select an atom in the four dimensions (4 estimated parameters).

Atresidual = opts.A'*residual;

% Choose the params that satify the following condition.
[~,idx] = max(abs(real(Atresidual)));
%idx = idx(1,1);

% the four variables's constrained region
% Ax <= b where x = [x1;x2;x3;x4]
%
A = [ -1  0 0 0;1 0 0 0 ; 0 -1 0 0;0 1 0 0;0 0 -1 0; 0 0 1 0; 0 0 0 -1; 0 0 0 1];
b = [-opts.B(1,1,1); opts.B(2,1,1) ;-opts.B(1,2,1); opts.B(2,2,1); -opts.B(1,1,2); opts.B(2,1,2) ;-opts.B(1,2,2); opts.B(2,2,2)];

% Minimize the " min_scal_prod " function in the constrained region.
% - param = [x1; x2; x3; x4]
%
fObj = @(param) min_scal_prod(residual,param,opts.atom,opts.datom);
options = optimoptions(@fmincon,'Display','off','GradObj','on','Algorithm','sqp');
param_new = fmincon(fObj,[opts.param_grid(1,idx);opts.param_grid(2,idx);opts.param_grid(3,idx); opts.param_grid(4,idx)],A,b,[],[],[],[],[],options);
val = scal_prod( residual , param_new , opts.atom );
out.p_screen = 0;

out.param_new = param_new;
out.val = val;

end
