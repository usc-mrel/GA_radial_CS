function x = CG_solver(x0,param,lambda,  a)

% x = argmin || NUFFT( sMaps * x ) - kU ||^2 + lambda*|| x - a ||^2

% initial image
n_iter = 0;
max_iter = 20;
x = x0;
r0 = param.E' * (param.y - param.E * x0)+lambda*(a-x0);
r = r0;
p = r0;
err = 1;

% if param.max_iter ~= 0
%     max_iter = param.max_iter;
% end


while (1)
    q = param.E' * (param.E *p)+lambda*p;
    alpha = r(:)' * r(:) / (p(:)' * q(:) );
    x = x + alpha * p;
    r_new = r - alpha * q;
    
    beta = r_new(:)' * r_new(:) / (r(:)' * r(:));
    p = r_new + beta * p;
    
    r = r_new;
    err_old = err;
    err= r(:)' * r(:) / (r0(:)' * r0(:));
    n_iter = n_iter + 1;
    

    fprintf('\n iteration #%d, alpha = %f, beta = %f, error = %f', n_iter, alpha, beta, abs(err));

    if (n_iter > max_iter) %|| ( abs(err) < tol*abs(err_old) )
        break;
    end
end
    

