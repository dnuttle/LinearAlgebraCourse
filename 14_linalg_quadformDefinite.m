%%
%     COURSE: Linear algebra: theory and implementation
%    SECTION: Quadratic form and definiteness
%      VIDEO: The quadratic form in algbebra
% Instructor: mikexcohen.com
% 
%%

% matrix and vector
S = [ 1 3 -2; 0 3 4; -5 -2 4 ];
w = [ -2 4 3 ]';

% compute the quadratic form
qf = w'*S*w;


n = length(w); % used for plotting

% show the matrices
figure(1), clf
subplot(131)
imagesc(S)
axis square, axis off
title('Matrix S')

subplot(132)
imagesc(w)
axis square, axis off
set(gca,'xlim',[.5 n+.5]-n/3)
title('Vector w')

subplot(133)
imagesc(qf)
axis square, axis off
set(gca,'xlim',[.5 n+.5]-n/3,'ylim',[.5 n+.5]-n/3)
title('Quadratic form: w^TSw')


%%
%     COURSE: Linear algebra: theory and implementation
%    SECTION: Quadratic form and definiteness
%      VIDEO: The quadratic form in geometry
% Instructor: mikexcohen.com
% 
%%

% some different matrices
S{1} = [ 4 4; 4 9 ];
S{2} = [ -4 -1; 3 -5 ];
S{3} = [ 0 1; 2 0 ];
S{4} = [ 1 1; 1 1 ];

% range for vector w
n = 30;
wRange = linspace(-2,2,n);

% initialize quadratic form matrix
qf = zeros(length(wRange),length(wRange));

figure(2), clf
for i=1:length(S)
    
    
    % compute QF
    for xi=1:n
        for yi=1:n
            
            % this w
            w = [ wRange(xi) wRange(yi) ]';
            
            % QF
            qf(xi,yi) = w'*S{i}*w;
        end
    end
    
    % show the surface
    subplot(2,2,i)
    surf(wRange,wRange,qf')
    
    % make it look a bit nicer
    shading interp, axis square
    set(gca,'fontsize',12,'clim',[-1 1]*max(abs(qf(:)))*.6)
    xlabel('w_1'), ylabel('w_2'), zlabel('Quad.form energy')
end


colormap jet


%%
%     COURSE: Linear algebra: theory and implementation
%    SECTION: Quadratic form and definiteness
%      VIDEO: Eigenvectors and the quadratic form surface
% Instructor: mikexcohen.com
% 
%%

% create a happy little symmetric matrix
A = [1 2; 2 3];

% weights along each dimension
xi = -2:.1:2;

quadform = zeros(length(xi));
for i=1:length(xi)
    for j=1:length(xi)
        
        % define the vector
        x = [xi(i) xi(j)]';
        
        % compute the quadratic form of matrix A
        quadform(i,j) = x'*A*x/(x'*x);
        
    end
end

figure(3), clf, hold on

% draw the surface
surf(xi,xi,quadform'), shading interp
title('Visual representation of quadratic form of matrix A')
xlabel('W_1'), ylabel('W_2'), zlabel('energy')
rotate3d on, axis square

% eigendecomposition
[V,D] = eig(A);
V = V*2; % stretch for visualization

% plot the eigenvectors
zfact = max(quadform(:));
plot3([0 V(1,1)],[0 V(2,1)],[1 1]*zfact,'m','linew',5)
plot3([0 V(1,2)],[0 V(2,2)],[1 1]*zfact,'k','linew',5)

colormap jet


%%
%     COURSE: Linear algebra: theory and implementation
%    SECTION: Quadratic form and definiteness
%      VIDEO: Application of the normalized quadratic form: PCA
% Instructor: mikexcohen.com
% 
%%


% simulation parameters
N = 1000; % time points
M =   20; % channels

% time vector (radian units)
t = linspace(0,6*pi,N);

% relationship across channels (imposing covariance)
chanrel = sin(linspace(0,2*pi,M))';

% create mean-centered dataset
data = bsxfun(@times,repmat( sin(t),M,1 ),chanrel) + randn(M,N)/3;
data = bsxfun(@minus,data,mean(data,2));


% compute covariance matrix
covmat = data*data'/(N-1);


% eigendecomposition of covariance matrix
[evecs,evals] = eig( covmat );


% sort eigenvalues and eigenvectors
[evecs,evals] = eig( covmat );
[evals,soidx] = sort( diag(evals),'descend' );
evecs = evecs(:,soidx);

% convert eigenvalues to percent variance explained
evals = 100*evals/sum(evals);


% compute component time series
r = 2; % two components
comp_time_series = evecs(:,1:r)'*data;


% visualize and interpret the results
figure(4), clf

% eigenvalues
subplot(231)
plot(evals,'s-','linew',2,'markerfacecolor','w')
axis square
xlabel('Component number'), ylabel('\lambda (% total variance)')

% eigenvectors
subplot(232)
plot(evecs(:,1:2),'s-','linew',2,'markerfacecolor','w')
axis square
xlabel('Channel'), ylabel('PC weight')
legend({'PC1';'PC2'})

% original channel modulator
subplot(233)
plot(chanrel,'s-','linew',2,'markerfacecolor','w')
axis square
xlabel('Channel'), ylabel('PC weight')

% component time series
subplot(212)
plot(1:N,comp_time_series)
xlabel('Time (a.u.)'), ylabel('Activity')
legend({'PC1';'PC2'})


%%
%     COURSE: Linear algebra: theory and implementation
%    SECTION: Quadratic form and definiteness
%      VIDEO: Matrix definiteness, geometry, and eigenvalues
% Instructor: mikexcohen.com
% 
%%


% some different matrices
S{1} = [ 4 4; 4 9 ];
S{2} = [ -4 -1; -3 -5 ];
S{3} = [ 0 1; 2 0 ];
S{4} = [ 1 1; 1 1 ];
S{5} = [ -1 -2; -3 -6 ];

% range for vector w
n = 30;
wRange = linspace(-2,2,n);

% initialize quadratic form matrix
qf = zeros(length(wRange),length(wRange));

figure(5), clf
for i=1:length(S)
    
    
    % compute QF
    for xi=1:n
        for yi=1:n
            
            % this w
            w = [ wRange(xi) wRange(yi) ]';
            
            % QF
            qf(xi,yi) = w'*S{i}*w;
        end
    end
    
    % show the surface
    subplot(2,3,i)
    surf(wRange,wRange,qf')
    
    % make it look a bit nicer
    shading interp, axis square
    set(gca,'fontsize',12,'clim',[-1 1]*max(abs(qf(:)))*.6)
    xlabel('w_1'), ylabel('w_2'), zlabel('Quad.form energy')
    
    
    %% compute the matrix's definiteness based on the eigenvalues
    
    % get eigenvalues
    evals = eig(S{i});
    
    % we care about their signs
    esign = sign(evals);
    
    % test for signs (note: this test is valid only for 2x2 matrices!)
    if sum(esign)==2
        defcat = 'Positive definite';
    elseif sum(esign)==1
        defcat = 'Positive semidefinite';
    elseif sum(esign)==0
        defcat = 'Indeterminant';
    elseif sum(esign)==-1
        defcat = 'Negative semidefinite';
    elseif sum(esign)==-2
        defcat = 'Negative definite';
    end
    
    % add title
    title([ defcat ' matrix' ])
    
end


colormap jet

%%

