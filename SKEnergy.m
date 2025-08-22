function E = SKEnergy(N,J,h)

Y = ff2n(N);
n = size(Y,1);
E = zeros(n,1);

for i=1:n
    y = Y(i,:);
    x = 2.*y - 1;
    
    E(i,1) = (x * J * x') + h*sum(x);
    
    % E = 0;
    % for i=1:N
    %     for j=1:N
    %         E = E + x(i)*x(j)*J(i,j);
    %     end
    % end
    % E = E/sqrt(N);
end

end