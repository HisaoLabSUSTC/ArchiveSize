function [index] = CDR(pop,n)
    [N, M] = size(pop);
    % the previous solution
    Pre = zeros(M,N); 
    % the next solution
    Next = zeros(M,N); 
    % crowding distance
    crowd = zeros(N,1); 
    % precalculate
    for i=1:M
        [~, I] = sort(pop(:,i));
        for j=1:N
            if j==1
                Pre(i,I(j))=inf;
                crowd(I(j)) = inf;
            else
                Pre(i,I(j))=I(j-1);
                crowd(I(j)) = crowd(I(j))+pop(I(j),i)-pop(I(j-1),i);
            end
            if j==N
                Next(i,I(j))=inf;
                crowd(I(j)) = inf;
            else
                Next(i,I(j))=I(j+1);
                crowd(I(j)) = crowd(I(j))+pop(I(j+1),i)-pop(I(j),i);
            end
        end
    end
    
    % remove (N-n) solutions based on the crowding distances
    ind_sel = ones(1,N);
    for k=1:N-n
        [~,I]=min(crowd);
        ind_sel(I)=0;
        crowd(I)=inf;

        for i=1:M
            % update the crowding distances of the previous solution
            % and the next solution.
            A = Pre(i,I);
            B = Next(i,I);
            if A~=inf
                if B==inf
                    crowd(A)=inf;
                else
                    crowd(A)=crowd(A) - pop(I,i) +  pop(B,i);
                end
            end
        
            if B~=inf
                if A==inf
                    crowd(B)=inf;
                else
                    crowd(B)=crowd(B) + pop(I,i) -pop(A,i);
                end
            end
            % update the link
            Next(i,A) = B;
            Pre(i,B) = A;
        end
    end
    
    index = 1:N;
    index(ind_sel==0)=[];
end

