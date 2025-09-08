function M=unfoldope(Atom,nb,wsoc)
    % Input:
    % Atom={'Atom1','Atom2',...}: atoms in order as same as in .win file. Ex: Atom={'Tb','Mn','O','O','O',....}
    % nb=[n1,n2,....]: number of orbits for each atom in order
    % wsoc : 1 for w/ soc, 0 for w/o soc
    % Note:
    % 1. require length(nb)==length(Atom)
    % 
    % Output:
    % M: Operator for unfolding spectral weight A=-imag(trace(M*G)), where G=1/(E-H)
    
    %check
    if length(Atom)~=length(nb)
        error('Require length(Atom)==length(nb).')
    end
    % label Atom
    n=0;
    b={};
    for iA=1:length(Atom)
        b=cat(2,b,{n+(1:nb(iA))});
        n=n+nb(iA);
    end
    % build M
    M=zeros(sum(nb));
    AA=unique(Atom);
    for iA=1:length(AA)
        mk=strcmp(AA{iA},Atom);
        nn=nb(mk);
        bb=cell2mat(b(mk));
        M(bb,bb)=kron(ones(sum(mk)),eye(nn(1)));
    end
    M=kron(eye(wsoc+1),M);
end