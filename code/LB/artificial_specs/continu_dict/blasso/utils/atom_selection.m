function [ out ] = atom_selection( residual , opts )

Atresidual = opts.A'*residual;
[~,idx] = max(abs(real(Atresidual)));
idx = idx(1,1);
fObj = @(param) min_scal_prod2(residual,param,opts.atom,opts.datom);
A = [ -1  0 ;1 0  ; 0 -1 ;0 1];
b = [-opts.B(1,1); opts.B(2,1) ;-opts.B(1,2); opts.B(2,2)];
%lb = [opts.B(1,1);opts.B(1,2)];
%ub = [opts.B(2,1); opts.B(2,2)];
options = optimoptions(@fmincon,'Display','off','GradObj','on','Algorithm','sqp');
param_new = fmincon(fObj,[opts.param_grid(1,idx);opts.param_grid(2,idx)],A,b,[],[],[],[],[],options);
val = scal_prod( residual , param_new , opts.atom );
out.p_screen = 0;

out.param_new = param_new;
out.val = val;

end
