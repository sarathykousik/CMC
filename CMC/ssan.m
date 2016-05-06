%%%% Modifed by Chaitanya

% 
% synth is the original signal to be filtered
% Z is the reconstructed signal
% r is the difference between recosntructed and original
% vr is the relative norm difference


function [y,r,vr]=ssan(x1,L)


% Step1 : Build trayectory matrix

   N=length(x1); 
   N
   if L>N/2;L=N-L;end
	K=N-L+1; 
   X=zeros(L,K);  
	for i=1:K
	  X(1:L,i)=x1(i:L+i-1); 
	end
    
% Step 2: SVD

   S=X*X'; 
	[U,autoval]=eig(S);
	[d,i]=sort(-diag(autoval));  
   d=-d;
   U=U(:,i);sev=sum(d); 
	plot((d./sev)*100),hold on,plot((d./sev)*100,'rx');
	title('Singular Spectrum');xlabel('Eigenvalue Number');ylabel('Eigenvalue (% Norm of trajectory matrix retained)')
   V=(X')*U; 
   rc=U*V';

% Step 3: Grouping

   I=input('Choose the  components to reconstruct the series in the form I=[i1,i2:ik,...,iL]  ')
   Vt=V';
   rca=U(:,I)*Vt(I,:);

% Step 4: Reconstruction

   y=zeros(N,1);  
   Lp=min(L,K);
   Kp=max(L,K);

   for k=0:Lp-2
     for m=1:k+1;
      y(k+1)=y(k+1)+(1/(k+1))*rca(m,k-m+2);
     end
   end

   for k=Lp-1:Kp-1
     for m=1:Lp;
      y(k+1)=y(k+1)+(1/(Lp))*rca(m,k-m+2);
     end
   end

   for k=Kp:N
      for m=k-Kp+2:N-Kp+1;
       y(k+1)=y(k+1)+(1/(N-k))*rca(m,k-m+2);
      end
   end

   figure;subplot(2,1,1);hold on;xlabel('Data poit');ylabel('Original and reconstructed series')
   plot(x1);grid on;plot(y,'r')
       
   
   r=x1-y';
    subplot(2,1,2);plot(r,'g');xlabel('Data poit');ylabel('Residual series');grid on
   vr=(sum(d(I))/sev)*100;
 