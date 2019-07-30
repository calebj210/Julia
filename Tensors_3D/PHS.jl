### Polyharmonic-spline defintions
ϕ(X,Y,m) = norm(X-Y)^m;
ϕ_xi(X,Y,m,i) = m*(X[i]-Y[i])*norm(X-Y)^(m-2);
ϕ_xii(X,Y,m,i) =
    begin
        nm = norm(X-Y);
        m*nm^(m-2)+m*(m-2)*(X[i]-Y[i])^2*nm^(m-4);
    end
ϕ_xij(X,Y,m,i,j) = m*(m-2)*(X[i]-Y[i])*(X[j]-Y[j])*norm(X-Y)^(m-4);
