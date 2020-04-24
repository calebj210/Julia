### Polyharmonic-spline defintions
ϕ(X,Y,m) = norm(X-Y)^m;
ϕ_xi(X,Y,m,i) = (m-2) >= 0 ? m*(X[i]-Y[i])*norm(X-Y)^(m-2) : 0.0;
ϕ_xii(X,Y,m,i) =
    begin
        nm = norm(X-Y);
        tmp = 0.0;
        if m-2 >= 0 
            tmp += nm^(m-2);
        end
        if m-4 >= 0
            tmp += (m-2)*(X[i]-Y[i])^2*nm^(m-4);
        end
        tmp *= m;
        
        return tmp
    end
ϕ_xij(X,Y,m,i,j) = m*(m-2)*(X[i]-Y[i])*(X[j]-Y[j])*norm(X-Y)^(m-4);
