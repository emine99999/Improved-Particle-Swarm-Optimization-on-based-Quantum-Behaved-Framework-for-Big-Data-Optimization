function y_obj = Objective(x,Amatrix,Mixed,ICAcomponent,DtypeG,SingMul,PopSize,coll)%#coder
y_obj= zeros(SingMul,PopSize);
for i=1:PopSize
    
    s1=reshape(x(i,:),coll,DtypeG)'; %% reshape x to be a mtrix of DtypeG*col
    X1 = Amatrix* s1;
    %% ==================== Pearson correlation coefficient ===============

    COR1 =corr_man2(X1',Mixed',coll);

    %%Minimize f1 = 1?(N × M) ?ij (Sij ? S1ij)2
    musum = sum(sum((ICAcomponent - s1).^2));
    %%Minimize f2 = 1?(N2 ? N) ?i,j ? i (Cij2) + 1?N ?i (1 ? Cii)2
    %% part 1
    f2_part1= (sum(sum(COR1.^2))-sum(diag(COR1.^2)))/(DtypeG^2-DtypeG);
    %% part2
    f_2part2= sum(diag((1-COR1).^2))/DtypeG;
    
    %% switch between single and multi-objective cases
    switch (SingMul)
        case 2
            y_obj(1,i) = f2_part1+f_2part2;
            y_obj(2,i) = musum/ (DtypeG*coll);
        case 1
            y_obj(1,i) = f2_part1+f_2part2+ (musum/ (DtypeG*coll)); %%f1+f2
    end
   
end
end

function COR1 = corr_man2 (A,B,n)

xc = bsxfun(@minus,A,sum(A,1)/n);  % Remove mean

    yc = bsxfun(@minus,B,sum(B,1)/n);  % Remove mean
    COR1 = xc' * yc; % 1/(n-1) doesn't matter, renormalizing anyway
    dx = sqrt(sum(abs(xc).^2, 1)); % sqrt first to avoid under/overflow
    dy = sqrt(sum(abs(yc).^2, 1)); % sqrt first to avoid under/overflow
    COR1 = bsxfun(@rdivide,COR1,dx'); COR1 = bsxfun(@rdivide,COR1,dy); % coef = coef ./ dx'*dy;


end