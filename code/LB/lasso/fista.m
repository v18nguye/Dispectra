function [ x , fc_lasso , fc_lassoDual ] = fista( y , opts)

x = opts.xinit;
if(opts.disp)
    disp('Fista running...')
end

t = 1/(2*opts.L);
alpha = opts.lambda*t;
tFISTA = 1;
xFISTA = x;
fc_lasso = zeros(1,opts.maxIter);
fc_lassoDual = zeros(1,opts.maxIter);

for iter = 1 : opts.maxIter
        
    xold = x;    
    x = xFISTA+t*opts.A'*(y-opts.A*xFISTA);
    x = lasso_thsld(x,alpha);
    
    fc_lasso(iter) = lasso_primal(y,opts.A,x,opts.lambda);
    
    dualv = (y-opts.A*x)/norm(opts.A'*(y-opts.A*x),'inf');
    fc_lassoDual(iter) = .5*norm(y)^2-.5*norm(y-opts.lambda*dualv,2)^2;
    dualgap = fc_lasso(iter)-fc_lassoDual(iter);
    if(dualgap<=opts.tol)
        fc_lasso = fc_lasso(1:iter-1);
        fc_lassoDual = fc_lassoDual(1:iter-1);
        break;
    end
    
    tFISTA_old = tFISTA;
    tFISTA = (1+sqrt(1+4*tFISTA^2))/2;
    xFISTA = x+(tFISTA_old-1)/tFISTA*(x-xold);
    
end

if(opts.disp)
    disp('============')
    if(iter<opts.maxIter)
        disp('Convergence reached')
    else
        disp('Maximum iterations reached')
    end
    disp(['Value of the lasso primal function: ',num2str(fc_lasso(end))])
    disp('============')
end

end


function lp = lasso_primal( y , A , x , lambda )
lp = .5*norm(y-A*x,2)^2+lambda*norm(x,1);
end


function thsld = lasso_thsld( x , alpha )
if(isreal(x))
    thsld = max(abs(x)-alpha,0).*sign(x);
else
    thsld = max(abs(x)-alpha,0).*exp(1i*angle(x));
end
end
