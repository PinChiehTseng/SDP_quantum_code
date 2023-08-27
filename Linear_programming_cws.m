warning('on','all')

value = [];
for n = 10
    for d = 4%2: n
        count = [];
        for K = 6: (2^(n))
            b = 0;
            clear i j k;
    
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
                    variable x(1, n+1);
                    variable y(1, n+1);
        
                    con = [];
        
                    ob = cvx(zeros(1));
                    ob_2 = cvx(zeros(1));
                    for i = 0: n
                        ob = ob + x(i+1);
                        ob_2 = ob_2 + y(i+1);
                    end
                    clear i;
        
                    maximize(ob)
                    subject to
                    x(1) == 1;
                    y(1) == 1;
                    ob_2 == K*(2^n);
                    for j = 0: n
                        x(j+1) >= 0;
                        y(j+1) >= 0;
                        x(j+1) <= y(j+1);
                    end
        
                    clear j;
                    
                    for p = 0: n
                        tot = cvx(zeros(1));
                        for l = 0: n
                            tot = tot + krav(l, n, p) * x(l+1);
                        end
                        if p < d
                            (K / (2^n))*tot == x(p+1);
                        end
                        (K / (2^n))*tot == y(p+1);
                        con = [con, tot];
                    end
        
                    clear p d_1 D tot;
        
                    for p = 0: n
                        tot = cvx(zeros(1));
                        for l = 0: n
                            tot = tot + ((-1)^l)*krav(l, n, p) * x(l+1);
                        end
                        %((K^2) / (2^n))*tot >= 0;
                        %con = [con, tot];
                    end
        
                    clear p d_1 D tot;
        
                cvx_end
                if cvx_status == "Infeasible" || cvx_status == "Overdetermined"
                   count = [count, K];
                   b = 1;
                end
                if K == (2^(n))
                   count = [count, (2^n)+1];
                end
                if b == 1
                    break;
                end
        end
        value = [value,"("+string(n)+","+string(min(count)-1)+","+string(d)+")"];
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