warning('on','all')

value = [];
sh = 1;
mi = 1;
c = 4;
Q = 2^c;
for n = c+1: 10
    for d = 3: n
        for K = 1: 2^(n+c)
                sB = K*(2^(n-c));
    
                cvx_clear
                cvx_begin sdp
                    cvx_precision('default')
                    %cvx_precision('high')
                    %cvx_precision('best')
                    %cvx_precision('medium')
                    %cvx_precision('low')
                    %cvx_solver sdpt3
                    %cvx_solver_settings( 'maxit', 1000 );
                    %cvx_solver SeDuMi
                    cvx_solver MOSEK
                    %cvx_solver_settings -clearall;
                    %variable x(1, n+1) integer;
                    %variable y(1, n+1) integer;
                    variable x(n+1, c+1);
    
        
                    ob = 0;
                    for i = 0: n
                        ob = ob + x(i+1, 1);
                    end
                    clear i;
        
                    maximize(ob)
                    subject to
                    x(1, 1) == 1;
                    for i = 0: n
                        for j = 0: c
                            x(i+1, j+1) >= 0;
                            if i == 0 && j > 0
                                x(i+1, j+1) == 0;
                            end
                        end
                    end
        
                    clear j i;

                    oby = 0;
                    obt = 0;
                    for p = 0: n
                        row = [];
                        for q = 0: c
                            tot = 0;
                            for l = 0: n
                                for r = 0: c
                                    tot = tot + (K / (2^(n+c)))*krav(l, n, p) * krav(r, c, q) * x(l+1, r+1);
                                end
                            end
                            if q == 0 && p < d
                                tot == x(p+1, 1);
                                if p == 0
                                    tot == 1;
                                end
                            end
                            tot >= x(p+1, q+1);
                            if q == 0
                                oby = oby + tot;
                            end
                            obt = obt + tot;
                        end
                    end
                    oby == sB;
                    obt == K*(2^(n+c));
                    clear p d_1 D tot q l r z tot_1 row;
                    if sh == 1 || mi == 1
                        enum = [];
                        for z = 0: n+c
                            tot = 0;
                            for p = 0: n
                                for q = 0: c
                                    if p+q == z
                                        tot = tot + x(p+1, q+1);
                                    end
                                end
                            end
                            enum = [enum, tot];
                        end
    
                        clear z p q tot;
                    end


                    if sh == 1
                        for p = 0: n+c
                            tot = 0;
                            for l = 0: n+c
                                tot = tot + ((-1)^l)*krav(l, n+c, p) * enum(l+1);
                            end
                            (K / (2^(n+c)))*tot >= 0;
                        end
                        clear p tot l;
                    end

                    if mi == 1
                        for p = 0: n+c
                            tot = 0;
                            for l = 0: n+c
                                tot = tot + (K / (2^(n+c)))*krav(l, n+c, p) * enum(l+1);
                            end
                            tot >= 0;
                        end
                        clear p tot l;
                    end
                cvx_end
                if cvx_status == "Infeasible"
                    A = "("+string(n)+","+string(K-1)+","+string(d)+")";
                    value = [value, A];
                    break
                end
                if K == 2^(n+c) && cvx_status == "Solved"
                    A = "("+string(n)+","+string(K)+","+string(d)+")";
                    value = [value, A];
                    break
                end
        end
    end
end

disp(value)



function A = comb(a, b)
    ans = 0;
    if a >= 0 && b >= 0 && b <= a
        ans = nchoosek(a, b);
    end
    A = ans;
end

function A = krav(x, n, k)
    ans = 0;
    for j = 0: k
        ans = ans + (((-1)^j))*(3^(k-j))*comb(x,j)*comb(n-x, k-j);
    end
    A = ans;
end