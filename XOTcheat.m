function XOTcheat()

% All primal and dual SDPs for the non-DI scenarios
% Takes ~10seconds to compute all results when using the simplification that some registers in the cheating-Bob scenario can be assumed to be classical

global basis pauli msqops msqpsi
%Convention: basis{1} = |0> = spin +Z, basis{2} = |1> = spin -Z
basis = {[1; 0], [0; 1]};
pauli = {[0 1; 1 0], [0 -i; i 0], [1 0; 0 -1]};
rng(1); 
addpath(genpath('mosek'))
addpath(genpath('YALMIP-master'))

%msqops is labelled such that (:,:,a,b,j) corresponds to outcome j of the measurement in row a, column b
msqops = zeros(2^2,2^2,3,3,2); %This version is defined from the Hermitian observables; code has also been verified using the projectors specified in the paper
msqops(:,:,1,1,1)=(eye(4)+Tensor({pauli{3},eye(2)}))/2;msqops(:,:,1,2,1)=(eye(4)+Tensor({eye(2),pauli{3}}))/2;msqops(:,:,1,3,1)=(eye(4)+Tensor({pauli{3},pauli{3}}))/2;
msqops(:,:,2,1,1)=(eye(4)+Tensor({eye(2),pauli{1}}))/2;msqops(:,:,2,2,1)=(eye(4)+Tensor({pauli{1},eye(2)}))/2;msqops(:,:,2,3,1)=(eye(4)+Tensor({pauli{1},pauli{1}}))/2;
msqops(:,:,3,1,1)=(eye(4)-Tensor({pauli{3},pauli{1}}))/2;msqops(:,:,3,2,1)=(eye(4)-Tensor({pauli{1},pauli{3}}))/2;msqops(:,:,3,3,1)=(eye(4)+Tensor({pauli{2},pauli{2}}))/2;
for a=1:3
for b=1:3
    if rank(msqops(:,:,a,b,1))~=2 || norm(msqops(:,:,a,b,1)*msqops(:,:,a,b,1) - msqops(:,:,a,b,1))>0
        error('Issue with magic-square projectors')
    else
        msqops(:,:,a,b,2) = eye(4) - msqops(:,:,a,b,1);
    end
end
end

%Sanity test: Should reproduce magic-square probabilities
msqpsi=(ket([0 0 0 0])+ket([0 1 0 1])+ket([1 0 1 0])+ket([1 1 1 1]))/2; %System ordering qX1 qX2 qY1 qY2
msqprobs=zeros(2^3,2^3,3,3);
for a=1:3
for b=1:3
for x=1:8
for y=1:8
    xvec=dec2bin(x-1,3)-'0';yvec=dec2bin(y-1,3)-'0'; %str to num witchcraft
    opA=msqops(:,:,a,1,xvec(1)+1)*msqops(:,:,a,2,xvec(2)+1)*msqops(:,:,a,3,xvec(3)+1);
    opB=msqops(:,:,1,b,yvec(1)+1)*msqops(:,:,2,b,yvec(2)+1)*msqops(:,:,3,b,yvec(3)+1);
    if norm(opA*opA-opA)>0 || norm(opB*opB-opB)>0
        error('Issue with game outcome projectors')
    end
    msqprobs(x,y,a,b)=msqpsi'*Tensor({opA,opB})*msqpsi;
%     [num2str([xvec yvec a b]) '  ' num2str([msqprobs(x,y,a,b)])] %Displays input/output tuples and corresponding probabilities
    if norm(msqprobs(x,y,a,b) - 1/8*(mod(sum(xvec),2)==0 && mod(sum(yvec),2)==1 && xvec(b)==yvec(a)))>0 %Some type abuse here by treating logical value as floating point
        error('Incorrect output distribution')
    end
end
end 
end
end

% Cheating Alice, untrusted state, primal
sigma=sdpvar(2^2,2^2,'hermitian','complex'); 
Qb=sdpvar(2^4,2^4,3,'hermitian','complex');
constr = [sigma>=0, trace(sigma)==1];

sumQb = zeros(2^4); fobj=0;
for b=1:3
    sumQb = sumQb + Qb(:,:,b);
    fobj = fobj + trace(choiB(b)*Qb(:,:,b)')/3;
    constr=[constr, Qb(:,:,b)>=0];
end
constr=[constr, sumQb==Tensor({eye(4),sigma})];
fobj = -fobj; %YALMIP runs as minimisation by default

fobj
constr
options=sdpsettings('verbose',1,'solver','mosek','savesolveroutput',1);
sol = optimize(constr,fobj,options);

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

sigmasol=value(sigma); Qbsol=value(Qb);
cheatAu=-value(fobj)
sigmasol
format short
if norm(imag(sigmasol))+norm(imag(Qbsol(:)))<10^-12 %Removes small imaginary components
    sigmasol=real(sigmasol);
    Qbsol=real(Qbsol);
end
smallinds=find(abs(Qbsol)<10^-12);Qbsol(smallinds)=0; %Removes small components
Qbsol

% Cheating Alice, untrusted state, dual
yalmip('clear')
lambda=sdpvar(1); R=sdpvar(2^4,2^4,'hermitian','complex');

fobj=lambda
constr=[lambda*eye(4)-PartialTrace(R,[1],[4 4])>=0 ,R>=choiB(1)/3, R>=choiB(2)/3, R>=choiB(3)/3]
sol = optimize(constr,fobj,sdpsettings('verbose',1,'solver','mosek','savesolveroutput',1));

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

cheatAudual=value(fobj)
Rsol=value(R);
format short
if norm(imag(Rsol(:)))<10^-12 %Removes small imaginary components
    Rsol=real(Rsol);
end
smallinds=find(abs(Rsol)<10^-12);Rsol(smallinds)=0; %Removes small components
Rsol

% Cheating Alice, trusted state, primal
yalmip('clear')
Eb=sdpvar(4*4,4*4,3,'hermitian','complex');
constr=[];

fobj=0;
sumEb=zeros(4*4);
for b=1:3
    rhoAb = zeros(16); eigstates=umeasB(b);
    for y=1:4
        cpart = zeros(4); cpart(y,y)=1; %Classical part of Alice's side-info, i.e. Bob's message
        eigstate=eigstates(y,:); proj=eigstate'*eigstate; 
        qpart = PartialTrace(Tensor({eye(4),proj})*msqpsi*msqpsi'*Tensor({eye(4),proj}),[2],[4 4]); %Quantum part of Alice's side-info, i.e. Bob's message. This is a subnormalised state that automatically incorporates the probability of outcome y.
        rhoAb = rhoAb + Tensor({cpart,qpart}); 
    end
    fobj = fobj + (1/3)*trace(Eb(:,:,b)*rhoAb); 
    constr = [constr, Eb(:,:,b) >= 0];
    sumEb = sumEb + Eb(:,:,b);
end
constr = [constr, sumEb==eye(4*4)];

fobj=-fobj; %YALMIP runs as minimisation by default
options=sdpsettings('verbose',1,'solver','mosek','savesolveroutput',1);
sol = optimize(constr,fobj,options);

format long
cheatAt = -value(fobj)
format short

% % Cheating Alice, trusted state, dual
yalmip('clear')
sigmaA = sdpvar(16,16,'hermitian','complex');
constr = [sigmaA >= 0];
%Constructing cq state storing Bob's measurement choice and Alice's side-info
for b=1:3
    rhoAb = zeros(16); eigstates=umeasB(b);
    for y=1:4
        cpart = zeros(4); cpart(y,y)=1; %Classical part of Alice's side-info, i.e. Bob's message
        eigstate=eigstates(y,:); proj=eigstate'*eigstate; 
        qpart = PartialTrace(Tensor({eye(4),proj})*msqpsi*msqpsi'*Tensor({eye(4),proj}),[2],[4 4]); %Quantum part of Alice's side-info, i.e. Bob's message. This is a subnormalised state that automatically incorporates the probability of outcome y.
        rhoAb = rhoAb + Tensor({cpart,qpart}); 
    end
    constr=[constr, (1/3)*rhoAb <= sigmaA];
end

fobj = trace(sigmaA)
constr
options=sdpsettings('verbose',1,'solver','mosek','savesolveroutput',1);
sol = optimize(constr,fobj,options);

format long
cheatAtdual = value(trace(sigmaA)) 
format short

% Cheating Bob, trusted/untrusted state (iteration 0/1 in the loop below), primal then dual
PcompYA1 = Tensor({ket([1 0])*ket([1 0])',diag([1 0 0])}) + Tensor({ket([0 1])*ket([0 1])',diag([0 1 0])}) + Tensor({ket([0 0])*ket([0 0])',diag([0 0 1])}) + + Tensor({ket([1 1])*ket([1 1])',eye(3)});
ErejXYA1G = Tensor({ket([0 0])*ket([0 0])' , Tensor({PcompYA1,eye(4)}) + Tensor({eye(4*3)-PcompYA1,eye(4)-ket([0 0])*ket([0 0])'})});
for x1=1:2
for x2=1:2
for x1p=1:2
for x2p=1:2
    if not(isequal([x1 x2],[1 1])) && not(isequal([x1 x2],[x1p x2p]))
        stateG=ket([x1p x2p]-ones(1,2)); stateX=ket([x1 x2]-ones(1,2));
        ErejXYA1G = ErejXYA1G + Tensor({stateX*stateX',eye(4*3),stateG*stateG'});
    end
end
end
end
end
if norm(ErejXYA1G*ErejXYA1G-ErejXYA1G)>10^-15
    error('Issue with Mrej')
end

cU = zeros(4*4*3*4); %System ordering qX Y A1 G
for a=1:3
    stateA=zeros(3,1); stateA(a)=1; 
    cU = cU + Tensor({umeasA(a),eye(4),stateA*stateA',eye(4)});
end

if norm(imag(cU))>0 || norm(imag(ErejXYA1G))>0 %Verifying these operators have no imaginary components, so it's safe to assume rho1,rho2 real in the SDP
    error('Imaginary components')
end

for untrusted=0:1
% Primal
yalmip('clear')
dims1 = [4 4]; rho1q = sdpvar(prod(dims1),prod(dims1),'symmetric','real'); %System ordering qX Y
dims2 = [4 4 3 4]; rho2q = sdpvar(prod(dims2),prod(dims2),'symmetric','real'); %System ordering qX Y A1 G
% Simplify matrices by projecting classical registers onto computational basis 
rho1=zeros(prod(dims1));
for basisnum=1:4
    proj=zeros(4);proj(basisnum,basisnum)=1;proj=Tensor({eye(4),proj});
    rho1=rho1+proj'*rho1q*proj;
end
rho2=zeros(prod(dims2));
for basisnum=1:prod(dims2(2:4))
    proj=zeros(prod(dims2(2:4)));proj(basisnum,basisnum)=1;proj=Tensor({eye(4),proj});
    rho2=rho2+proj'*rho2q*proj;
end
constr = [rho1>=0, rho2>=0];

fobj = trace(cU'*(eye(4*4*3*4)-ErejXYA1G)*cU*rho2);
fobj = -fobj; %YALMIP runs as minimisation by default

constr = [constr, PartialTrace(rho2,[4],dims2)==Tensor({rho1,eye(3)/3})];
if untrusted == 0
    constr = [constr, PartialTrace(rho1,[2],dims1)==eye(4)/4]; %Trusted-state version
else
    constr = [constr, trace(rho1)==1]; %Untrusted-state version 
end

fobj
constr
options=sdpsettings('verbose',1,'solver','mosek','savesolveroutput',1);
sol = optimize(constr,fobj,options);

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

if untrusted == 0
    cheatBt=-value(fobj)
else
    cheatBu=-value(fobj)
end

rho1sol=value(rho1);
smallinds=find(abs(rho1sol)<10^-12);rho1sol(smallinds)=0;
rho2sol=value(rho2);
smallinds=find(abs(rho2sol)<10^-12);rho2sol(smallinds)=0;
% format short
% rho1sol
% rho2sol
% dlmwrite('rho1sol.csv',rho1sol,'precision','%.15f')
% dlmwrite('rho2sol.csv',rho2sol,'precision','%.15f')

% Dual
yalmip('clear')
if untrusted == 0
    dims1 = [4 4 3]; Q1 = sdpvar(prod(dims1),prod(dims1),'symmetric','real'); %System ordering qX Y A1
    dims2 = [4]; Q2 = sdpvar(prod(dims2),prod(dims2),'symmetric','real'); %System ordering qX
    
    fobj=trace(Q2/4)
    constr=[1/3*PartialTrace(Q1,[3],dims1)<=Tensor({Q2,eye(4)}), cU'*(eye(4*4*3*4)-ErejXYA1G)*cU<=Tensor({Q1,eye(4)})]
    options=sdpsettings('verbose',1,'solver','mosek','savesolveroutput',1);
    sol = optimize(constr,fobj,options);
    
    cheatBtdual=value(fobj)
else
    eta = sdpvar(1);
    dims = [4 4 3]; Q = sdpvar(prod(dims),prod(dims),'symmetric','real'); %System ordering qX Y A1

    fobj = eta
    constr=[1/3*PartialTrace(Q,[3],dims)<=eta*eye(4*4), cU'*(eye(4*4*3*4)-ErejXYA1G)*cU<=Tensor({Q,eye(4)})]
    options=sdpsettings('verbose',1,'solver','mosek','savesolveroutput',1);
    sol = optimize(constr,fobj,options);
    
    cheatBudual=value(fobj)
end

end %End untrusted/trusted loop

summary = [cheatAu cheatAudual cheatAt cheatAtdual cheatBt cheatBtdual cheatBu cheatBudual]

end
%______________________________________________________________________________

function choiB = choiB(b)
% Choi matrix (unnormalized), following the convention that the map is applied to the first system
global msqops
choiB = zeros(2^4);
for y1=1:2
for y2=1:2
    measop=msqops(:,:,1,b,y1)*msqops(:,:,2,b,y2);
    if norm(measop*measop-measop)>0 || min(eigs(measop))<-10^-15 || rank(measop)~=1
        error('Issue with Bob projectors')
    end
    [eigvec eigval] = eigs(measop,1);
    choiB = choiB + Tensor({ket([y1 y2]-ones(1,2))*ket([y1 y2]-ones(1,2))',eigvec*eigvec'});
end
end
end

function umeasA = umeasA(a)
%Each row corresponds to an eigenvector of the PVM
global msqops
umeasA = zeros(2^2);
for x1=1:2
for x2=1:2
    measop=msqops(:,:,a,1,x1)*msqops(:,:,a,2,x2);
    if norm(measop*measop-measop)>0 || min(eigs(measop))<-10^-15 || rank(measop)~=1
        error('Issue with Alice projectors')
    end
    [eigvec eigval] = eigs(measop,1);
    umeasA = umeasA + ket([x1 x2]-ones(1,2))*eigvec';
end
end
end

function umeasB = umeasB(b)
%Each row corresponds to an eigenvector of the PVM
global msqops
umeasB = zeros(2^2);
for y1=1:2
for y2=1:2
    measop=msqops(:,:,1,b,y1)*msqops(:,:,2,b,y2);
    if norm(measop*measop-measop)>0 || min(eigs(measop))<-10^-15 || rank(measop)~=1
        error('Issue with Bob projectors')
    end
    [eigvec eigval] = eigs(measop,1);
    umeasB = umeasB + ket([y1 y2]-ones(1,2))*eigvec';
end
end
end

function ket = ket(bitlist)
global basis
qubits = cell(length(bitlist),1);
for n = 1:length(bitlist)
    qubits{n} = basis{bitlist(n) + 1};
end
ket = Tensor(qubits);
end