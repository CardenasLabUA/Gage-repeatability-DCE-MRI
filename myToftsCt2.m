function [c_toi] = myToftsCt2(ktrans,kep,t,Cp)
% Simulate concentration using the standard Tofts model
% [c_toi] =myToftsCt2(ktrans,kep,t,Cp)
% ktrans=pars(1);
% kep=pars(2);
% t=X(:,1);
% Cp=X(:,2);

n_points=length(t);
% expo=zeros(1,n_points);
% crpexp=expo;
% c_toi=crpexp;

for k = 2:n_points
    int_t = t(k);

    for j = 1:k
        dummy_t = t(j);
        expo(j) =exp(-((kep.*(int_t-dummy_t))));
        crpexp(j) = Cp(j)*expo(j);
    end
    t2 = t(1:k);
    %%
    crpexp_integral = trapz(t2,crpexp);
    %     crpexp_integral = trapzfm(t2,crpexp); % quicker?
    c_toi(k) = ktrans*crpexp_integral;
end
