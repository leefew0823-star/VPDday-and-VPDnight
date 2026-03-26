function [cof, r2, pvl, k] = myRR_vif(Y, X, vif_limit, k_step, k_limit)
[n, p] = size(X);

mx = mean(X);
my = mean(Y);
stdx = std(X);
stdy = std(Y);

stdx(stdx == 0) = 1e-8;
stdy(stdy == 0) = 1e-8;

Z = (X - mx(ones(n,1),:)) ./ stdx(ones(n,1),:);
Zy = (Y - my) ./ stdy;

R = corrcoef(Z);

kz = 0;
while kz <= k_limit

    inv_mat = inv(R + kz * eye(p));
    vif_ridge = diag(inv_mat * R * inv_mat);

    if max(vif_ridge) < vif_limit
        break;
    end
    kz = kz + k_step;
end
k = kz;

cof = inv(R + k * eye(p)) * (Z' * Zy / (n - 1));


Zy_pred = Z * cof;


SS_res = sum((Zy - Zy_pred).^2);
SS_tot = sum((Zy - mean(Zy)).^2);


r2 = 1 - (SS_res / SS_tot);


df_model = p;
df_error = n - p - 1;


if SS_res < 1e-10
    F_stat = inf;
    pvl = 0;
else
    F_stat = ((SS_tot - SS_res) / df_model) / (SS_res / df_error);
    pvl = 1 - fcdf(F_stat, df_model, df_error);
end

end