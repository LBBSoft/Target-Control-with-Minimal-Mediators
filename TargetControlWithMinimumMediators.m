clear all
clc
tic

% Generate scale-free network with parameters (N,mo,m)
N=100;
mo=5;
m=5;
hwait = waitbar(0,'Please wait. Generating directed scale-free adjacency matrix');
A = zeros(N);
E = 0;
for i=1:mo
	for j=i+1:mo
		
		A(i,j) =1;
		A(j,i) = 1;
		E = E + 2;
	end
    %% 
end
% Second add remaining nodes with a preferential attachment bias - rich get
% richer
for i=mo+1:N
	waitbar(i/N,hwait,sprintf('Please wait. Generating directed scale-free adjacency matrix\n%.1f %%',(i/N)*100));
	curr_deg =0;
	while(curr_deg<m)
		sample = setdiff(1:N,[i find(A(i,:))]);
		j = datasample(sample,1);
		b = sum(A(j,:))/E;
		r = rand(1);
		if(b>r)
			r = rand(1);
			if(b>r)
			A(i,j) = 1;
			A(j,i) = 1;
			E = E +2;
			else
			A(i,j) = 1;
			E = E +1;
			end
		else
			no_connection = 1;
			while(no_connection)
				sample = setdiff(1:N,[i find(A(i,:))]);
				h = datasample(sample,1);
				b = sum(A(h,:))/E;
				r = rand(1);
				if(b>r)
					r = rand(1);
					if(b>r)
					A(h,i) = 1;
					A(i,h) = 1;
					E = E +2;
					else
					A(i,h) = 1;
					E = E+1;
					end
					no_connection = 0;
				end
			end
		end
		curr_deg = sum(A(i,:));
	end
end

% A is state matrix
GG=A;
A=A';
n=length (A(:,1));
CC=(1:1:n);

% C is target nodes

% % full control
% C=(1:1:n);
% s=length (C(:,1));

% % Random selection scheme
% I=(1:1:n)';
% Percent=10;
% c=false(size(I));
% d=fix(Percent/100*numel(I));
% index=randi(numel(I),[1 d]);
% c(index)=true;
% h=nnz(c);
% while h<d
%     index=randi(numel(I),[1 d-h]);
%     c(index)=true;
%     h=nnz(c);
% end
% 
% C=I(c);
% s=length (C(:,1));

% Local selection scheme
I=(1:1:n)';
Percent=10;
d=fix(Percent/100*numel(I));
e=1;
while e<d
    bb=fix(n*rand(1))+1;
    B=[bb];
    z=[bb];
    s=[bb];
    while (~isempty(s))
        b=z(length(z));
        adj=find(A(b,:)==1);
        for t=1:length(adj)
            if ~ismember(adj(t),B)
                B=[B adj(t)];
            end
        end
        s=setdiff(B,z);
        if ~isempty(s)
            z=[z s(1)];
        end
    end
    e=size(z);
end
C=z(1:d);
s=length (C(:,1));

% construct set of accessible node with different lengths for each vertex

W=cell(n);
for i=1:n
    W{i,1}=i;
    for j=2:n;
        f=length(W{i,j-1});
        for t=1:f
            x=W{i,j-1}(t);
            for k=1:n
                if A(k,x)==1;
                    if ismember (k, W{i,j})
                    else
                        W{i,j}=[W{i,j},k];
                    end
                end
            end
        end
    end
end

% construct intersect of accessible node set with target node

for i=1:n
    for j=1:n
        H{i,j}=intersect(W{i,j},C);
    end
end

%  construct matrix that determines each vertex with what length walk from any vertex is available

K=cell(n);
for i=1:n;
    for t=1:n;
        for j=1:n;
            if ismember (t, H{i,j})
                K{t,i}=[K{t,i},j];
            end
        end
    end
end

% construct aggregation target nodes that are accessible for each vertex

Z=cell(n,1);
for i=1:n;
    for j=1:n;
        if isempty( H{i,j});
        else
            PPP=setdiff( H{i,j},Z{i});
            uy=length(H{i,j});
            PPPP=[];
            for ff=1:uy;
                if ismember(H{i,j}(ff) , PPP)
                    PPPP=[PPPP,H{i,j}(ff)];
                end
            end
            Z{i}=[Z{i},PPPP];
        end
    end
end

% construct upper bound for controlled set of each vertex

Y=cell(n,1);
for i=1:n;
    if isempty(Z{i});
    else
        p=length(Z{i});
        Y{i}(1)=Z{i}(1);
        q=length(Y{i});
        for t=2:p;
            for j=1:q;
                if   isempty ( setxor( K{Z{i}(t),i}, K{Y{i}(j),i}));
                    break
                elseif j==q
                    Y{i}=[Y{i},Z{i}(t)];
                    q=length(Y{i});
                end
            end
        end
    end
end

% construct controllable nodes set of each vetex

for i=1:n
    l=length(Y{i});
    G=zeros(n,n);
    ll=l;
    for j=1:n;
        for t=1:n;
            if  ismember(t,H{i,j})
                G(t,j)=1;
            end
        end
    end
    vv=0;
    KK=setdiff(CC,Y{i});
    G(KK,:)=[];
    u=rank(G);
    v=l-u;
    vv=vv+v;
    while v~=0;
        l=u;
        G=G([1:u],:);
        u=rank(G);
        v=l-u;
        vv=vv+v;
    end
    Y{i}=Y{i}(1:ll-vv);
end
YY=Y;

% determining driver nodes
for i=1:n
    M(i)=length(Y{i});
end
[a,b]= max(M);
if ismember(b ,setdiff(CC,C));
    for i=1:n;
        if ismember(i,C) && M(i)==M(b);
            b=i;
        end
    end
end
driver=[];
O=[];
driver=[driver,b];
X=cell(n,2);
O=[O,Y{b}];
X{1,1}=b;
X{1,2}=Y{b};
tt=2;
N=C;
N=setdiff(N,O);
while length(N)>=1;
    for i=1:n;
        if i~=b;
            UU{i}=setdiff(Y{i},Y{b});
            mmm=length(Y{i});
            U=[];
            for u=1:mmm
                if ismember (Y{i}(u), UU{i});
                    U=[U,Y{i}(u)];
                end
            end
            Y{i}=U;
        end
    end
    Y{b}=[];
    for i=1:n;
        M(i)=length(Y{i});
    end
    [a,b]= max(M);
    if ismember(b ,setdiff(CC,C));
        for i=1:n;
            if ismember(i,C) && M(i)==M(b);
                b=i;
            end
        end
    end
    X{tt,1}=b;
    X{tt,2}=Y{b};
    tt=tt+1;
    driver=[driver,b];
    O=[O,Y{b}];
    N=setdiff(N,O);
end
mm=length(driver);
XZ=[];
for i=1:mm;
    ll=mm-i+1;
    L=[];
    for j=1:mm
        if j~=ll ;
            if ismember(j, XZ);
            else
                L=[L, YY{X{j,1}}];
            end
        end
    end
    if isempty(setdiff(X{ll,2} , L));
        XZ=[XZ,ll];
        k=length(X{ll,2});
        for t=1:k;
            for j=1:mm;
                if j~=ll;
                    if ismember(j, XZ);
                    else
                        if ismember (X{ll,2}(t) , YY{X{j,1}})
                            X{j,2}=[X{j,2}, X{ll,2}(t)];
                            break
                        end
                    end
                end
            end
        end
    end
end
driver([XZ])=[];
X([XZ],:)=[];

% determining the ratio of the number of the mediator nodes to the number of the target nodes 
mm=length(driver);
mm;
P=cell(mm,n);
for i=1:mm
    P{i,1}=X{i,1};
    m=length(X{i,2});
    for j=2:m+1
        P{i,j}=X{i,2}(j-1);
    end
    
end 
for i=1:mm
    P{i,1}=[P{i,1};0];
    m=length(X{i,2});
    s=1;
    for j=2:m+1
        for t=s:n
            if ismember(P{i,j} , H{X{i,1},t});
                P{i,j}=[P{i,j};t];
                break
            end
        end
        s=t;
    end
end
sd=0;
for i=1:mm;
    m=length(X{i,2});
    sd=sd+P{i,2}(2);
    for j=3:m+1;
        for ll=2:j-1;
            l=j+1-ll;
            if ismember((P{i,j}(2)-P{i,l}(2)+1) , K{P{i,j}(1) , P{i,l}(1)} );
                sd=sd + (P{i,j}(2)-P{i,l}(2)) ;
                break
            end 
        end
    end
end
cc=length (C(1,:));
sdc=sd/cc;

driver
RatioOfMediatorToTarget=sdc

toc

