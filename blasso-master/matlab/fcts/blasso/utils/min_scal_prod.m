function [ msp , gradsp ] = min_scal_prod( residual , param , atom , datom )
msp = -sign(real(atom(param)'*residual))*real(atom(param)'*residual);
gradsp = -sign(real(atom(param)'*residual))*real(datom(param)'*residual);
end
