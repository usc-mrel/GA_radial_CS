function [samp] = Goldenratiosampling(nx,nt,nrays);

ang=0;
 for frameno = 1:nt,
    
    %  x=-nx/2:1:nx/2-0.5;  % Create an array of N points between -N/2 and N/2
      x=linspace(-0.5,0.5,nx); 
     i=1;
     while i<=nrays
         klocn=complex(x*cos(ang),x*sin(ang));
         kloc_all(:,i)=klocn;
         %ang = ang + (sqrt(5)-1)/2*pi; % the golden ratio step
        ang = ang+2*pi/(3.23606797); %
   if ang > 2*pi
      ang = mod(ang,2*pi); % always be between 0 and 2pi
   end
    i= i +1;
   end
    
    samp(:,:,frameno)=kloc_all;
 end
 
%samp = 0.5*(samp./max(samp(:))); 