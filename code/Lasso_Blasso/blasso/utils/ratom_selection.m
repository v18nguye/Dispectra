function [ out ] = ratom_selection( residual , opts )

Atresidual = opts.A'*residual;
[~,idx] = max(abs(real(Atresidual)));
idx = idx(1,1);
fObj = @(param) r_min_scal_prod2(residual,param,opts.atom,opts.datom);
%C = [ -1 ; 1 ];
%b = [ -opts.B(1); opts.B(2) ];
A = [ -1  0 0 0;1 0 0 0 ; 0 -1 0 0;0 1 0 0;0 0 -1 0; 0 0 1 0; 0 0 0 -1; 0 0 0 1];
b = [-opts.B(1,1,1); opts.B(2,1,1) ;-opts.B(1,2,1); opts.B(2,2,1); -opts.B(1,1,2); opts.B(2,1,2) ;-opts.B(1,2,2); opts.B(2,2,2)];
%lb = [opts.B(1,1);opts.B(1,2)];
%ub = [opts.B(2,1); opts.B(2,2)];
options = optimoptions(@fmincon,'Display','off','GradObj','on','Algorithm','sqp');
param_new = fmincon(fObj,[opts.param_grid(1,idx);opts.param_grid(2,idx);opts.param_grid(3,idx); opts.param_grid(4,idx)],A,b,[],[],[],[],[],options);
val = scal_prod( residual , param_new , opts.atom );
out.p_screen = 0;

out.param_new = param_new;
out.val = val;

end