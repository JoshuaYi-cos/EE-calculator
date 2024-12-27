function [entropy, dos] = d1_rand_timelike(L, points, tstep, fill)
    % Finite 1d chain, calculate the entanglement entropy of the subsystem
    % with support over the point 0. 

    % fill the bottom energies by ordering the eigenvalues and finding the
    % correct permutation matrices to act on the vectors
    [basis, energies] = d1_H_rand_set(L);
    [energies_s,ind] = sort(diag(energies));

    basis_s = basis(:,ind); 
    
    % initialize the points (column vector)
    t = tstep*(0:1:(points-1))';
    
    % calculate B matrix
    % multiply the first row of V^\dagger with the first column of V such
    % that we grab the magnitude of the first component of every
    % eigenvector. First component since our region has support over 0.
    
    % Auxiliary matrices
    M = exp(1i.*energies_s*t');
    J = ones(1,points);

    % find the filling
    N = round(fill*L);

    % delta cutoff for eliminating small evalues
    delta = 1e-12;

    aux = zeros(1,L);

    for i = 1:L
        % calculate change of basis
        P = abs(basis_s(i,:)).^2';
        B = M'*((P*J).*M);

        % extract the M matrix from the diagonalized B matrix
        [V,D]=eig(B); 
        
        D = diag(D);
        indices = (D<delta);
        D = D(~indices);
        V = V(:,~indices);
        D = diag(D);
        
        % test out whether our diagonalization is okay, we could do O'*O-B
        O=diag(1./sqrt(diag(D)))*V';

        C = M(1:N,:)'*((P(1:N)*J).*M(1:N,:));

        d = conj(O)*C*O.';
        [Vt,Dt] = eig(d);
        Dt = diag(Dt);
        Dt = real(Dt); % get rid of complex numbers
        Dt = Dt(Dt~=0); % get rid of zeros
        Dt = Dt(Dt~=1); % get rid of ones
        S = -sum(Dt.*log(Dt) + (1-Dt).*log(1-Dt));
        aux(i) = real(S);
    end

    dos = basis_s;
    entropy = mean(aux);

    % old code

%     for i = 1:points
%         for j = i:points
%             B(i,j) = basis_s(1,:)*diag(exp(1i.*energies_s.*(t(i)-t(j))))*basis_s(1,:)';
%             B(j,i) = conj(B(i,j));
%         end
%     end
   
    
    
    % calculate correlation function
%     cor = zeros(points);
%     
%     for i = 1:points
%         for j = i:points
%             cor(i,j) = basis_s(1,1:N)*diag(exp(1i.*energies_s(1:N).*(t(i)-t(j))))*basis_s(1,1:N)';
%             cor(j,i) = conj(cor(i,j));
%         end
%     end
    
    % calculate the entropy
end