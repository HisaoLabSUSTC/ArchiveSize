function [B, ind_sel] = DSS(a,sel_num)
start = tic;
if size(a,1)<=sel_num
    B=a; 
    ind_sel = 1: size(a,1);
else
    %The data is already normalized
    anorm=a;
    [sz_a, M]=size(a);  
    [~, ind_extreme] = max(a(:,1));
    first_extreme=a(ind_extreme,:); %get the first extreme solution


    %initialize B as empty set 
    B=[]; 

    %initialize v and d 
    v=zeros(sz_a, 1);     %zeroes represent false 
    d=ones(sz_a,1)*inf;   %initiliaze distance to a very large number, i.e., inf

    %select first extreme solutions 
    B=[B;first_extreme]; 
    v(ind_extreme)=1; 

    for l=1:sz_a
        if v(l)==0
            d(l)=min(sqrt(sum((anorm(l,:)-anorm(ind_extreme,:)).^2)),d(l));
        end
    end

    sz_B=size(B,1); 
    ind_sel=[];
    ind_sel=[ind_sel;ind_extreme]; 
    %select an isolated solution to B 
    while sz_B < sel_num
        z=find(v==0); 
        max_distance=max(d(z)); 
        j=find(d==max_distance);
        j=j(1);
        B=[B;a(j,:)];
        v(j)=1; 
        d(j)=inf;

        for l=1:sz_a
            if v(l)==0
                d(l)=min(sqrt(sum((anorm(l,:)-anorm(j,:)).^2)),d(l));
            end
        end
        if toc(start)>3600
            break;
        end

        sz_B=size(B,1); 
        ind_sel=[ind_sel;j];
        %ind_sel=find(v==1); 
    end   
end
end

