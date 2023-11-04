warning('on','all')

clear all
sh = 1;
mi = 1;
slp = 1;
mat = 1;
sp_sh = 1;
values = [];
c = 8;

table = [
    4681, 1164, 307, 96, 24, 6, 3, 0];


for n = c+1: 9
    for d = 4%3: n
        count = 0;
        for K = table(n-c, d-2)
            if K == 0
                break
            end
            b = 0;
            order = 10^(-3);
            lst_1 = cell(n+1, 1);
            lst_2 = cell(n+1, 1);
            lst = [];
            sA = 2^(n-c);
            sB = K*(2^(n-c));
            for i = 1: n+1
                lst_1{i, 1} = ["0"];
                lst_2{i, 1} = [0, 0, 0, 0];
            end
            clear i;
            for i = 0: n
                for j = 0: n
                    for t = 0: min(i, j)
                        for p = 0: t
                            if per(i, j, t, p, n) ~= "N"
                                if ~ismember(per(i, j, t, p, n), lst_1{t-p+1, 1})
                                    k = t-p+1;
                                    lst_1{k, 1} = [lst_1{k, 1}, per(i, j, t, p, n)];
                                    lst_2{k, 1} = [lst_2{k, 1}; [i, j, t, p]];
                                    lst = [lst; [i, j, t, p]];
                                end
                            end
                        end
                    end
                end
            end
            clear i j t p;
            
            S = size(lst);
            num = S(1);
            disp(num*4 + (n+1)*(c+1));
            clear S;
            
            cvx_clear
            cvx_begin sdp
                cvx_precision('default')
                %cvx_precision('high')
                %cvx_precision('best')
                %cvx_precision('medium')
                %cvx_precision('low')
                %cvx_solver sdpt3
                %cvx_solver SeDuMi
                cvx_solver MOSEK
                variable x(num);
                variable y(num);
                variable y_1(num);
                variable y_2(num);
                if sh == 1 || mi == 1 || slp == 1
                    variable L(n+1, c+1);
                end
            
                object = cvx(zeros(1, n+1));
                for i = 0: n
                    A_in = find(i, 0, 0, 0, n, lst_1, lst_2);
                    B = check_num(A_in, lst);
                    object(1, i+1) = (3^i)*comb(n, i)*x(B);
                end
                clear i A_in B;
            
            
                maximize(sum(object))
                subject to
                %%definition condition
                A = find(0, 0, 0, 0, n, lst_1, lst_2);
                B = check_num(A, lst);
                x(B) == 1;
                y(B) == 1;
                y_1(B) == 1;
                y_2(B) == 1;
                clear A B;
                if mat == 1
                    for a = 0: n
                        for k = 0: floor(n/2)
                            if 2*k < n
                                disp([a, k])
                                [matrix_1, matrix_2] = generate(x, n, a, k+a, lst_1, lst_2, lst);
                                order*(matrix_2 - matrix_1) >= 0;
                                matrix_1*order >= 0;
                                [matrix_3, matrix_4] = generate(y, n, a, k+a, lst_1, lst_2, lst);
                                matrix_3*order >= 0;
                                order*(matrix_4 - matrix_3) >= 0;
                                [matrix_5, matrix_6] = generate(y_1, n, a, k+a, lst_1, lst_2, lst);
                                matrix_5*order >= 0;
                                order*(matrix_6 - matrix_5) >= 0;
                                [matrix_7, matrix_8] = generate(y_2, n, a, k+a, lst_1, lst_2, lst);
                                matrix_7*order >= 0;
                                order*(matrix_8 - matrix_7) >= 0;
                            end
                        end
                    end
                end
                clear a k;
                ob = 0;
                ob_y = [];
                ob_y_1 = 0;
                ob_y_2 = 0;
                for i = 0: n
                    A_In = find(i, 0, 0, 0, n, lst_1, lst_2);
                    B = check_num(A_In, lst);
                    ob = ob + (3^i)*comb(n, i)*x(B);
                    ob_y = [ob_y, (3^i)*comb(n, i)*y(B)];
                    ob_y_2 = ob_y_2 + (3^i)*comb(n, i)*y_2(B);
                    ob_y_1 = ob_y_1 + (3^i)*comb(n, i)*y_1(B);
                end
                clear i A_In B;
            
                ob*order == sB*order;
                sum(ob_y)*order == sA*order;
                ob_y_2*order == sA*order;
                ob_y_1*order == sB*order;
                %ob == ob_y_1;

                if slp == 1
                    L(1, 1) == 1;
                    for i = 0: n
                        for j = 0: c
                            L(i+1, j+1) >= 0;
                            if i == 0 && j > 0
                                L(1, j+1) == 0;
                            end
                        end
                    end
        
                    clear j i;
                    s_1 = 0;
                    s_2 = 0;
                    for p = 0: n
                        for q = 0: c
                            tot = 0;
                            tot_sh = 0;
                            for l = 0: n
                                for r = 0: c
                                    tot = tot + (K / (2^(n+c))) * krav(l, n, p) * krav(r, c, q) * L(l+1, r+1);
                                    tot_sh = tot_sh + (K / (2^(n+c))) * ((-1)^(l+r)) * krav(l, n, p) * krav(r, c, q) * L(l+1, r+1);
                                end
                            end
                            if sp_sh == 1
                                tot_sh+order >= 0;
                            end
                            if q == 0 && p < d
                                tot*order == order*L(p+1, 1);
                                if p == 0
                                    tot == 1;
                                end
                            end
                            tot*order >= order*L(p+1, q+1);
                            if q == 0 && p > 0
                                tot*order == object(1, p+1)*order;
                            end
                            if q == 0
                                s_1 = s_1 + tot;
                            end
                            s_2 = s_2 + tot;
                        end
                    end
                    order*s_1 == order*sB;
                    order*s_2 == order*K*(2^(n+c));
                    clear p d_1 D tot q l r z tot_sh;
                end

                if sh == 1 || mi == 1
                    enum = [];
                    for z = 0: n+c
                        tot = 0;
                        for p = 0: n
                            for q = 0: c
                                if p+q == z
                                    tot = tot + L(p+1, q+1);
                                end
                            end
                        end
                        enum = [enum, tot];
                    end
    
                    clear z p q tot;
                end

                if sh == 1
                    for p = 0: n+c
                        tot = cvx(zeros(1));
                        for l = 0: n+c
                            tot = tot + ((-1)^l)*krav(l, n+c, p) * enum(1, l+1);
                        end
                        order*tot >= 0;
                    end
                    clear p tot l;
                end

                if mi == 1
                    for p = 0: n+c
                        tot = cvx(zeros(1));
                        for l = 0: n+c
                            tot = tot + (K / (2^(n+c)))*krav(l, n+c, p) * enum(1, l+1);
                        end
                        order*tot >= 0;
                    end
                    clear p tot l;
                end
                if mat == 1
                for i = 0: n
                    disp(i);
                        for j = 0: n
                            for t = 0: min(i, j)
                                for p = 0: t
                                    if i+j-t-p >= 0 && i+j-t-p <= n
                                        A = find(i, j, t, p, n, lst_1, lst_2);
                                        B_1 = check_num(A, lst);
                                        B = find(i, 0, 0, 0, n, lst_1, lst_2);
                                        B_2 = check_num(B, lst);
                                        x(B_1) >= 0; 
                                        y(B_1) >= 0;
                                        y_1(B_1) >= 0;
                                        y_2(B_1) >= 0;
                                        -x(B_1) + x(B_2) >= 0;
                                        -y(B_1) + y(B_2) >= 0;
                                        -y_1(B_1) + y_1(B_2) >= 0;
                                        -y_2(B_1) + y_2(B_2) >= 0;
                                        y(B_1) <= y_1(B_1);
                                        y(B_1) <= y_2(B_1);
                                        y_1(B_1) <= x(B_1)*sB / sA;
                                        y_2(B_1) <= x(B_1)*sB / sA;
                                        y(B_1) <= x(B_1);
                                        if i < d && j < d
                                            y(B_1) == y_1(B_1);
                                            y_1(B_1) == y_2(B_1);
                                        end
                                        if i < d || j < d
                                            y_2(B_1) == y(B_1);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                clear i j t p B_1 B_2 A B;
            
            
            cvx_end

            if cvx_status == "Infeasible" || cvx_status == "Overdetermined"
               count = K;
            end
            if count ~= 0
                values = [values,"("+string(n)+","+string(count-1)+","+string(d)+")"];
                break
            end
        end
    end
end

disp(values)



function A = krav(x, n, k)
    ans = 0;
    for j = 0: k
        ans = ans + (((-1)^j))*(3^(k-j))*comb(x,j)*comb(n-x, k-j);
    end
    A = ans;
end

function B = check_num(A, lst)
    S = size(lst);
    num = S(1);
    for k = 1: num
        if lst(k, 1) == A(1) && lst(k, 2) == A(2) && lst(k, 3) == A(3) && lst(k, 4) == A(4)
            B = k;
        end
    end
end

function A = find(i, j, t, p, n, lst_1, lst_2)
    if per(i, j, t, p, n) ~= "N"
        D = size(lst_1{t-p+1, 1});
        for q = 1: D(2)
            C = lst_1{t-p+1, 1};
            if per(i, j, t, p, n) == C(q)
                B = lst_2{t-p+1, 1};
                A = B(q, :);
            end
        end
    end
end

function Y = per(i, j, t, p, n)
    O = "1";
    if i-t < 0 || j-t < 0 || i+j-t-p < 0 || i+j-t-p > n || i-p < 0 || j-p < 0 || t-p < 0
        O = "N";
    end
    if i > n || j > n
        O = "N";
    end
    V = perms([i j i+j-t-p]);
    a = [];
    item = 0;
    for q = 1: 6
        vector = V(q, :);
        m = vector(1, 1);
        l = vector(1, 2);
        k = vector(1, 3);
        check_1 = 0: n;
        if (ismember(k, check_1) && (ismember(m, check_1) && ismember(l, check_1)))
            if item > 0
                check_10 = 0;
                for F = 1: item
                    if a(F) == [m, l, k]
                        check_10 = 1;
                    end
                end
                if check_10 == 0
                    item = item + 1;
                    a = [a; [m, l, k]];
                end
            else
                item = item + 1;
                a = [m, l, k];
            end
        end
    end
    if item == 1
        A = a(1);
        B = a(2);
        C = a(3);
    end
    first = [];
    for q = 1: item
        first = [first, a(q, 1)];
    end
    one = min(first);
    f = [];
    f_num = 0;
    second = [];
    for r = 1: item
        if a(r, 1) == one
            f = [f; a(r, :)];
            f_num = f_num + 1;
            second = [second, a(r, 2)];
        end
    end
    if f_num == 1
        A = f(1);
        B = f(2);
        C = f(3);
    end
    two = min(second);
    g = [];
    g_num = 0;
    third = [];
    for s = 1: f_num
        if f(s, 2) == two
            g_num = g_num + 1;
            g = [g; f(s, :)];
            third = [third, f(s, 3)];
        end
    end
    if g_num == 1
        A = g(1);
        B = g(2);
        C = g(3);
    end
    three = min(third);
    h = [];
    h_num = 0;
    for l = 1: g_num
        if f(l, 3) == three
            h = f(l, :);
            h_num = 1;
        end
    end
    if h_num == 1
        A = h(1);
        B = h(2);
        C = h(3);
    end
    if O == "1"
        O = '['+string(A)+',' +string(B) +','+string(C) +']';
    end
    Y = O;
end

function B = comb_num(i, j, k, m, t)
    ans = 0;
    for u = 0:m
        ans = ans + ((-1)^(t-u))*comb(u,t)*comb(m-2*k,m-u-k)*comb(m-k-u,i-u)*comb(m-k-u,j-u);
    end
    B = ans;
end

function A = comb(a, b)
    ans = 0;
    if a >= 0 && b >= 0 && b <= a
        ans = nchoosek(a, b);
    end
    A = ans;
end

function A = alpha(i, j, t, p, a, k, n)
    B = 0;
    for g = 0: p
        z = a-g;
        B = B + ((-1)^z)*comb(a, g)*comb(t-a, p-g)*(2^(t-a-p+g));
    end
    A = comb_num(i-a, j-a, k-a, n-a, t-a)*(3^(((i+j)/2)-t))*B;
end

function A = enumerator(x, n, K, j, lst_1, lst_2, lst)
    B = 0;
    for i = 0: n
        C = find(i, 0, 0, 0, n, lst_1, lst_2);
        D = check_num(C, lst);
        B = B + (3^i)*comb(n, i)*krav(i, n, j)*x(D);
    end
    B = B / ((2^n)*K);
    A = B;
end

function [matrix_1, matrix_2] = generate(x, n, a, k, lst_1, lst_2, lst)
    goal_1 = cvx(zeros(n+a-2*k+1, n+a-2*k+1));
    goal_2 = cvx(zeros(n+a-2*k+1, n+a-2*k+1));
    for i = k: n+a-k
        for j = k: n+a-k
            a_1 = 0;
            b_1 = 0;
            for t = 0: min(i, j)
                for p = 0: t
                    if i+j-t-p >= 0 && i+j-t-p <= n
                        B = find(i, j, t, p, n, lst_1, lst_2);
                        A = check_num(B, lst);
                        a_1 = a_1 + alpha(i, j, t, p, a, k, n)*x(A);  
                        if i+j-t-p >= 0 && i+j-t-p <= n
                            B_1 = find(i+j-t-p, 0, 0, 0, n, lst_1, lst_2);
                            A_1 = check_num(B_1, lst);
                            b_1 = b_1 + alpha(i, j, t, p, a, k, n)*x(A_1);
                        end
                    end
                end
            end
            goal_1(i-k+1, j-k+1) = a_1;
            goal_2(i-k+1, j-k+1) = b_1;
        end
    end
    matrix_1 = goal_1;
    matrix_2 = goal_2;
 end