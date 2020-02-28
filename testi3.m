A=zeros(3,2,2);
A(1,:,:)=[1 0;
          0 1];
A(2,:,:)=[0 1;
          1 0];
A(3,:,:)=[1 -1;
          1  1]/sqrt(2);      

B=[1;1];

squeeze(A(3,:,:))*B

%uusi
% 0.7239         0    0.6289    0.1990
% 0.1068    0.4942    0.1060         0
% 0.1461    0.4952    0.1450         0

%vanha
% 0.7239         0    0.5984    1.7625
% 0.1122    3.8507    0.1060         0
% 0.1543    4.0335    0.1450         0


bifurcations=[1 2 3; 
              3 4 5];
%parameter contains parameters on all segments. First column contains Beta,
%second contains A0, 3. rho, 4. length of the segment,5. and 6. the
%starting and ending points of the segment (in cm). 7. is pext
% parameters=zeros(5,7);
% parameters(:,1)=[237.66; 303.796; 325.67; 399.714; 325.318];
% parameters(:,2)=[0.51; 0.106; 0.145; 0.031; 0.133];
% parameters(:,3)=rho;
% parameters(:,4)=[42.4; 23.5; 6.7; 7.9; 17.1];%lengths
% parameters(:,7)=pext;
% 
% seg=bifurcations(2,:)
% parameters(seg,2)
% A=parameters(seg,1)
% A(2)







