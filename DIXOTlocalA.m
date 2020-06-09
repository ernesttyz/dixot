function DIXOTlocalA()

%Estimated runtime 4.5hrs to compute all points

global sizeA sizeB sizeX
rng(1); 
tic;

npalocalA=1
npalocalB=1
sizeA=4;sizeB=4;sizeX=7;sizeY=3; %We will order Alice's measurements as [Mtest Mprot]

%opslist will be the list of operators indexing the rows/columns of the NPA matrix. Since all projectors are hermitian, we can use the same opslist to index both rows and columns.
%Each term opslist{j} represents some product of terms of the form P^(party)_{c|z}. 
%It is encoded as a matrix of the form [c1,c1,...; z1,z2,...] representing P_{c1|z1}P_{c2|z2}..., where z-values larger than sizeX index Bob's operators.
%Special: [] indicates the identity operator, which is consistent with using concatenation as operator multiplication

level1listA={}; level1listB={}; %Constructing lists of first-order operators
maxA=sizeA;maxB=sizeB; %Here we will use all the projectors, without any normalization-based eliminations
for x = 1:sizeX
    maxAfewer=sizeA; %Some of Alice's measurements have fewer outcomes
    if x>=4 
        maxAfewer=3;
    end
    for a = 1:maxAfewer
        level1listA{end+1} = [a;x];
    end
end
for y = 1:sizeY
    for b = 1:maxB
        level1listB{end+1} = [b;sizeX+y];
    end
end
%Constructing list of all products of order up to npalocalA for Alice; analogous list for Bob
for party=1:2
    if party==1
        level1list=level1listA; npalocal=npalocalA;
    else
        level1list=level1listB; npalocal=npalocalB;
    end
    opslist=[{ [] },level1list]; %Start with the identity operator and the first-order moments
    leveljlist = level1list; 
    for j=2:npalocal
        leveljnew={};
        for k=1:length(leveljlist)
            for l=1:length(level1list)
                temp=optype([leveljlist{k},level1list{l}]);
                if isequal(temp,0) || isequal(temp,[])
                    % Do nothing
                else
                    leveljnew{end+1}=temp;
                end
            end
        end
        leveljlist=leveljnew; opslist=[opslist,leveljlist];
    end
    if party==1
        opslistA=celluniques(opslist);
    else
        opslistB=celluniques(opslist);
    end
    clear('opslist','npalocal')
end %party

% Generates opslist as all products of opslistA x opslistB
opslist={ [] }; %Start with identity operator
for j=1:length(opslistA)
for k=1:length(opslistB)
    temp=optype([opslistA{j},opslistB{k}]);
    if isequal(temp,0) || isequal(temp,[])
        % Do nothing
    else
        opslist{end+1}=temp;
    end
end
end

opslist=celluniques(opslist);
numops=length(opslist)

% Performs the following reductions:
% (1) Moment type: Identity, zero, or other
% (2) Squared projectors are collapsed
% These reductions commute with each other since collapsing squared projectors does not change result of (1)
momentsref=cell(numops);
for j=1:numops
    for k=1:numops
        % This ordering assumes the operators from the row index j are left-multiplied to those from the col index k
        % Also, it assumes the row indices are the ones that are conjugated (in this context, that means the order of the operators is reversed first)
        momentsref{j,k} = optype([fliplr(opslist{j}),opslist{k}]);
    end
end

yalmip('clear')
momentsmat = sdpvar(numops,numops,'symmetric') %WLOG can assume all moments are real (as long as all constraints are real), since the complex conjugate of the gamma matrix yields another feasible point with the same values wrt objective and constraint.
constr=[];
posconstrs=0; normconstrs=0; nsconstrs=0;

datetime('now')
uniquemoms={ [] };uniquepos=[1,1]; %Assumes the top-left element is always the identity
if momentsref{1,1} ~= []
    momentsref{1,1}
    error('Top-left element is not identity!')
end
for j=1:numops
for k=1:numops 
    if isequal(momentsref{j,k},0)
        momentsmat(j,k)=0; momentsmat(k,j)=0;
    elseif isequal(momentsref{j,k}, [] )
        momentsmat(j,k)=1; momentsmat(k,j)=1;
    else
        duplicatepos=cell1Dpos(uniquemoms,momentsref{j,k});
        if length(duplicatepos)==0 %Have not seen this moment before
            uniquemoms{end+1}=momentsref{j,k};uniquepos=[uniquepos;j,k];
            
%             %Positivity constraints; apparently have negligible effect in this setting
%             op=momentsref{j,k}; posA=find(op(2,:)<=sizeX); posB=find(op(2,:)>sizeX);
%             opsA=op(:,posA); opsB=op(:,posB);
%             if isequal(fliplr(opsA),opsA) && isequal(fliplr(opsB),opsB) 
%                 constr=[constr,momentsmat(j,k)>=0];
%                 posconstrs=posconstrs+1;
%             end
            
        elseif length(duplicatepos)==1 %Have previously seen this moment
            rownum=uniquepos(duplicatepos(1),1); colnum=uniquepos(duplicatepos(1),2);
            momentsmat(j,k)=momentsmat(rownum,colnum); momentsmat(k,j)=momentsmat(rownum,colnum);
        else
            error('Issue with identifying duplicates!')
        end
    end
end
if mod(j,30)==0 %Prints timestamp every 30 rows
    j
    datetime('now')
end
end
constr=[constr,momentsmat>=0]; %NPA matrix is PSD

numuniques=length(uniquemoms) %Number of unique entries in the moment matrix
momentsmat
datetime('now')

% %Normalisation constraints; apparently have negligible effect in this setting
% for x=1:sizeX
% for y=1:sizeY
%     normterm = 0;
%     maxAfewer=sizeA;
%     if x>=4 
%         maxAfewer=3;
%     end
%     for a=1:maxAfewer 
%     for b=1:sizeB
%         pos=cell1Dpos(uniquemoms,[a b; x y+sizeX]);rownum=uniquepos(pos,1);colnum=uniquepos(pos,2);
%         if length(pos)==0
%             error('Element not found!')
%         end
%         normterm = normterm + momentsmat(rownum,colnum);
%     end
%     end
%     constr = [constr, normterm==1];
%     normconstrs=normconstrs+1;
% end
% end
% 
% %No-signalling constraints; apparently have negligible effect in this setting
% for x=1:sizeX
% maxAfewer=sizeA;
% if x>=4 
%     maxAfewer=3;
% end
% for a=1:maxAfewer 
%     pos=cell1Dpos(uniquemoms,[a; x]);rownum=uniquepos(pos,1);colnum=uniquepos(pos,2);
%     marginal=momentsmat(rownum,colnum);
% for y=1:sizeY
%     nsterm=0;
%     for b=1:sizeB
%         pos=cell1Dpos(uniquemoms,[a b; x y+sizeX]);rownum=uniquepos(pos,1);colnum=uniquepos(pos,2);
%         if length(pos)==0
%             error('Element not found!')
%         end
%         nsterm = nsterm + momentsmat(rownum,colnum);
%     end
%     constr = [constr, nsterm==marginal];
%     nsconstrs=nsconstrs+1;
% end
% end
% end
% for y=1:sizeY
% for b=1:sizeB
%     pos=cell1Dpos(uniquemoms,[b; y+sizeX]);rownum=uniquepos(pos,1);colnum=uniquepos(pos,2);
%     marginal=momentsmat(rownum,colnum);
% for x=1:sizeX
%     nsterm=0;
%     maxAfewer=sizeA;
%     if x>=4 
%         maxAfewer=3;
%     end
%     for a=1:maxAfewer 
%         pos=cell1Dpos(uniquemoms,[a b; x y+sizeX]);rownum=uniquepos(pos,1);colnum=uniquepos(pos,2);
%         if length(pos)==0
%             error('Element not found!')
%         end
%         nsterm = nsterm + momentsmat(rownum,colnum);
%     end
%     constr = [constr, nsterm==marginal];
%     nsconstrs=nsconstrs+1;
% end
% end
% end

testterm = {{ [] },[0]};
for x=1:3
for y=1:3
for a=1:4
for b=1:4
    avec=dec2bin(a-1,2)-'0';bvec=dec2bin(b-1,2)-'0';
    avec=[avec mod(avec(1)+avec(2),2)];
    bvec=[bvec mod(bvec(1)+bvec(2)+1,2)];
    if avec(y)==bvec(x)
%         [x y avec bvec]
        testterm = lincombplus(testterm, { {[a b; x y+sizeX]} , [1/9]});
    end
end
end
end
end

cheatterm = {{ [] },[0]};
for y=1:3
for b=1:4
    cheatterm = lincombplus(cheatterm, {{ [y b; b+3 y+sizeX] }, [1/3]} );
end
end

testvars=0;
ops=testterm{1};coeffs=testterm{2}
for j=1:length(coeffs)
    pos=cell1Dpos(uniquemoms,ops{j});rownum=uniquepos(pos,1);colnum=uniquepos(pos,2);
    if length(pos)==0
        ops{j}
        error('Element not found!')
    end
    testvars = testvars + coeffs(j)*momentsmat(rownum,colnum);
end
cheatvars=0;
ops=cheatterm{1};coeffs=cheatterm{2}
for j=1:length(coeffs)
    pos=cell1Dpos(uniquemoms,ops{j});rownum=uniquepos(pos,1);colnum=uniquepos(pos,2);
    if length(pos)==0
        ops{j}
        error('Element not found!')
    end
    cheatvars = cheatvars + coeffs(j)*momentsmat(rownum,colnum);
end

testproblist = (1:9)/10; 
resultslist = zeros(2,length(testproblist));
for pt = 1:length(testproblist)

testprob = testproblist(pt)
fobj = testprob*testvars + (1-testprob)*cheatvars

fobj = -fobj; %YALMIP runs as minimisation by default
% testprob
fobj=clean(fobj,10^-10)
% constr
length(constr)
posconstrs
normconstrs
nsconstrs
options=sdpsettings('verbose',1,'solver','mosek','savesolveroutput',1,'saveduals',0);
datetime('now')
sol = optimize(constr,fobj,options);
datetime('now')
cheatA=-value(fobj)
toc;

format long
if isfield(sol,'solveroutput')
    pobjval = sol.solveroutput.res.sol.itr.pobjval
    dobjval = sol.solveroutput.res.sol.itr.dobjval
    solstatus = sol.solveroutput.res.sol.itr.solsta
    absgap = pobjval-dobjval
    relgap = (pobjval-dobjval)/((pobjval+dobjval)/2)
else
    pobjval = []; dobjval = []; solstatus = [];
end

probssol=zeros(sizeA,sizeB,sizeX,sizeY);
for x=1:sizeX
for y=1:sizeY
for a=1:sizeA
for b=1:sizeB
    pos=cell1Dpos(uniquemoms,[a b; x y+sizeX]);rownum=uniquepos(pos,1);colnum=uniquepos(pos,2);
    if length(pos)==0
        probssol(a,b,x,y) = 0;
    else
        probssol(a,b,x,y) = value(momentsmat(rownum,colnum));
    end
end
end
end
end
probssol
smallinds = find(abs(probssol)<10^-5);
probssol(smallinds) = 0
cheatA=-value(fobj)

momentsmatvals = value(momentsmat);
save(strcat(num2str(round(100*testprob)),'_sols.mat'),'momentsmatvals')
resultslist(:,pt) = [testprob; cheatA]

toc;

end %end pt loop

end
%______________________________________________________________________________

function cell1Dpos = cell1Dpos(celllist,target)
cell1Dpos=[];
for j=1:length(celllist)
    if isequal(celllist{j},target)
        cell1Dpos= [cell1Dpos,j];
    end
end
end

function celluniques=celluniques(cellarray)
flattened=cellarray(:); celluniques={};
for j=1:length(flattened)
    isnew=1;
    for k=1:length(celluniques)
        if isequal(flattened{j},celluniques{k})
            isnew=0;
            break
        end
    end
    if isnew==1
        celluniques{end+1}=flattened{j};
    end
end
end

function reduceprojs = reduceprojs(lincomb)
% In each term of lincomb, eliminates all projectors corresponding to the "last" outcome of Alice or Bob
lincombnew2=reduce1stprojs(lincomb);
for j=1:100
    lincombold=lincombnew2; lincombnew2=reduce1stprojs(lincombnew2);
    if isequal(lincombnew2,lincombold)
        break
    elseif j==100
        lincombnew2
        lincombold
        j
        error('Reduction did not terminate!')
    end
end
reduceprojs=lincombnew2;
end

function reduce1stprojs = reduce1stprojs(lincomb)
% In each term of lincomb, eliminates the first projector corresponding to the "last" outcome of Alice or Bob
global sizeA sizeB sizeX
ops=lincomb{1};coeffs=lincomb{2};
if length(ops)~=length(coeffs)
    error('Length mismatch!')
end
for j=1:length(coeffs)
    opmatrix=ops{j};
    if length(opmatrix)==0 %Special case for the identity operator
        lincombnew = {{ [] }, [coeffs(j)]};
    else
        for k=1:size(opmatrix,2)
            if opmatrix(1,k)==sizeA && opmatrix(2,k)<=sizeX %Checks if projector is the "last" one of Alice
                opsnew=cell(1,sizeA); coeffsnew=coeffs(j)*[-ones(1,sizeA-1) 1];
                if length(opmatrix)==1
                    opsnew{sizeA}=[];
                else
                    temp=opmatrix; temp(:,k)=[]; opsnew{sizeA}=temp;
                end
                for a=1:sizeA-1
                    temp=opmatrix; temp(1,k)=a; opsnew{a}=temp;
                end
                lincombnew = {opsnew, coeffsnew};
                break
            elseif opmatrix(1,k)==sizeB && opmatrix(2,k)>sizeX %Checks if projector is the "last" one of Bob
                opsnew=cell(1,sizeB); coeffsnew=coeffs(j)*[-ones(1,sizeB-1) 1];
                if length(opmatrix)==1
                    opsnew{sizeB}=[];
                else
                    temp=opmatrix; temp(:,k)=[]; opsnew{sizeB}=temp;
                end
                for b=1:sizeB-1
                    temp=opmatrix; temp(1,k)=b; opsnew{b}=temp;
                end
                lincombnew = {opsnew, coeffsnew};
                break
            elseif k==size(opmatrix,2) %reached end of loop without finding a last-outcome projector
                lincombnew = {{ opmatrix }, [coeffs(j)]};
            end
        end %k loop
    end
    if j==1
        lincombfinal=lincombnew;
    else
        lincombfinal=lincombplus(lincombfinal,lincombnew);
    end
end %j loop
reduce1stprojs=lincombfinal;
end

function lincombtimes = lincombtimes(lincomb1,lincomb2)
ops1=lincomb1{1};coeffs1=lincomb1{2};ops2=lincomb2{1};coeffs2=lincomb2{2};
if length(ops1)~=length(coeffs1) || length(ops2)~=length(coeffs2)
    error('Length mismatch!')
end
opsnew={};coeffsnew=[];
for j=1:length(ops1)
for k=1:length(ops2)
    op=optype([ops1{j},ops2{k}]);
    if isequal(op,0)~=1
        opsnew{end+1}=op;
        coeffsnew(end+1)=coeffs1(j)*coeffs2(k);
    end
end
end
lincombtimes=lincombplus({opsnew,coeffsnew/2},{opsnew,coeffsnew/2}); %Lazy way to strip out zeros and duplicates
end

function lincombplus = lincombplus(lincomb1,lincomb2)
ops1=lincomb1{1};coeffs1=lincomb1{2};ops2=lincomb2{1};coeffs2=lincomb2{2};
%%%%%%%%%%%%%%%%%% The part between these lines makes lincombplus always merge projectors. Comment this out if that is not desired.
tempops=cell(1,2); bothops={ops1,ops2};
for term=1:2
    for opnum=1:length(bothops{term})
        tempops{term}{end+1}=mergeprojs(bothops{term}{opnum});
    end
end
ops1=tempops{1};ops2=tempops{2};
%%%%%%%%%%%%%%%%%%
if length(ops1)~=length(coeffs1) || length(ops2)~=length(coeffs2)
    error('Length mismatch!')
end
opsnew = celluniques([ops1 ops2]); coeffsnew=zeros(1,length(opsnew)); nonzeros=[];
for j=1:length(opsnew)
    %In the following, poslist may be empty (i.e. opsnew{j} does not appear in ops1/ops2), but it still works out since at least one of vals1,vals2 will be nonempty.
    poslist=cell1Dpos(ops1,opsnew{j});vals1=coeffs1(poslist); 
    poslist=cell1Dpos(ops2,opsnew{j});vals2=coeffs2(poslist);
    coeffsnew(j) = sum([vals1 vals2]);
    if coeffsnew(j)~=0
        nonzeros = [nonzeros,j];
    end
end
opsnew = opsnew(nonzeros); coeffsnew=coeffsnew(nonzeros);
if length(nonzeros)==0 %It might not be critical to check this, but I'd like to avoid singularities with lincomb={{},[]} if possible, by essentially defining a "canonical" zero element in the next line
    lincombplus = {{ [] },[0]}
else
    lincombplus = {opsnew, coeffsnew};
end
end

function optype = optype(opmatrix)
% Returns 0 if it contains consecutive terms of orthogonal projectors, returns mergeprojs(opmatrix) otherwise
temp=mergeprojs(opmatrix); 
for j=1:size(temp,2)-1
    temp1=temp(:,j);temp2=temp(:,j+1);
    if temp1(2)==temp2(2) && temp1(1)~=temp2(1) 
        temp=0;
        break
    end
end
optype=temp;
end

function mergeprojs = mergeprojs(opmatrix)
% Gathers Alice and Bob's projectors, then detects squared projector terms and collapses them
global sizeX
if length(opmatrix)==0
mergeprojs=opmatrix;
else
inputindices=opmatrix(2,:);
posA=(inputindices<=sizeX);
posB=(inputindices>sizeX);
collapsed = [opmatrix(:,posA), opmatrix(:,posB)]; nosquares=0; %0 indicates there are squares, 1 indicates there are no squares
for loopcounter=1:100 %We'll bound it at maximum 100 iterations
    nosquares=1; %Temporarily set its value to 1, will be reset to 0 if a square is found
    for j=size(collapsed,2)-1:-1:1 %Warning: using length(collapsed') causes problems when there is only one projector in opmatrix, because in that case opmatrix is only 2x1 and length() gives the number of rows rather than columns.
        if isequal(collapsed(:,j),collapsed(:,j+1)) 
            collapsed(:,j)=[];
            nosquares=0;
        end
    end
    if nosquares==1
        break
    end
    if j==100
        error('Merging projectors failed!')
    end
end
mergeprojs = collapsed;
end
end