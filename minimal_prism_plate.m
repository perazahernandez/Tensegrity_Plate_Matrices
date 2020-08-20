function [N,CB,CS] = minimal_prism_plate(q,h,d)
% This code provides the node locations and connectivity matrices
% for a tensegrity plate formed by minimal regular 3-bar prisms
%
% Inputs:
% q: Complexity of the plate
% h: Thickness of the plate
% d: Diameter of the plate
%
% Outputs:
% N: Matrix whose columns are the position vectors of the nodes
% CB: Connectivity matrix of the bars
% CS: Connectivity matrix of the strings
%
% Notes:
% The connectivity matrices satisfy the relations B = N*CB' and S = N*CS'
% where B is a matrix whose columns are the vectors along the length of the
% bars and S is a matrix whose columns are the vectors along the length of
% the strings
%
% This code is part of the Supplemental Material of
% Shuhui Jiang, Robert E. Skelton, and Edwin A. Peraza Hernandez, 2020, 
% "Analytical equations for the connectivity matrices and node positions 
% of minimal and extended tensegrity plates,"
% International Journal of Space Structures, 
% Vol. 35(3), pp. 47-68. https://doi.org/10.1177/0956059920902375

%% 1.1 cover area and radius of prism unit
Area=3*(d/2)^2*cos(pi/6);% Area - cover area of plate
r0=sqrt(2*Area/(3*sqrt(3)*(9*(2-sqrt(3))*(q-1)^2+3*(q-1)+1)));% r0 - radius of prism unit


%% 1.2. connectivity matrix
% CB - bar connectivity matrix
% CR - gravity center connectivity matrix of bars
% CS - string connectivity matrix
% CSt - top string connectivity matrix
% CSb - bottom string connectivity matrix
% CSv - vertical string connectivity matrix
[CB,CS,CR,CSt,CSb,CSv]=connectivitymatrix(q);



%% 1.3 initial node location of plate
[N]=nodematrix(q,h,r0);% N - node matrix



%% 1.4 plot
[rows,cols]=find(CS'); % rows - node No. for each string
[rowb,colb]=find(CB'); % rowb - node No. for each bar
[rowst,colst]=find(CSt');[rowsb,colsb]=find(CSb');[rowsv,colsv]=find(CSv');  % rowst,rowsb,rowsv - node No. for top, bottom, vertical string
nb=size(CB,1);ns=size(CS,1); % nb, ns - bar number and string number
nst=size(CSt,1);nsb=size(CSb,1);nsv=size(CSv,1); % nst, nsb, nsv - number of top, bottom, vertical strings
figure;
% plot bars(black lines)
for j=1:nb
plot3(N(1,rowb(2*j-1:2*j)),N(2,rowb(2*j-1:2*j)),N(3,rowb(2*j-1:2*j)),'linewidth',2,'color','k','marker','.'); hold on; 
end
% plot top strings(red lines)
for j=1:nst
plot3(N(1,rowst(2*j-1:2*j)),N(2,rowst(2*j-1:2*j)),N(3,rowst(2*j-1:2*j)),'linewidth',1,'color','r','marker','.'); hold on;
end
% plot bottom strings(blue lines)
for j=1:nsb
plot3(N(1,rowsb(2*j-1:2*j)),N(2,rowsb(2*j-1:2*j)),N(3,rowsb(2*j-1:2*j)),'linewidth',1,'color','b','marker','.'); hold on;
end
% plot vertical strings(green lines)
for j=1:nsv
plot3(N(1,rowsv(2*j-1:2*j)),N(2,rowsv(2*j-1:2*j)),N(3,rowsv(2*j-1:2*j)),'linewidth',1.5,'color','g'); hold on;
end
xlabel('x'); ylabel('y'); zlabel('z');
grid off;
axis equal;
axis on;
set(gca,'fontsize', 15,'linewidth',1.15)
set(gca,'ticklength',1.2*get(gca,'ticklength'))
drawnow;

end

%*****************************
% sub function
%*********************************
%*******************************
% 1. connectivity matrix
%********************************
function [CB,CS,CR,CSt,CSb,CSv]=connectivitymatrix(q) 
nu=1+3*q*(q-1);%number of units
no=(q*6-3)*2;% nodes only connect with strings
nb=3*nu;%number of bars
% 1) for bars
CBT=sparse([-eye(nb);eye(nb);zeros(no,nb)]);CRT=sparse([.5*eye(nb);.5*eye(nb);zeros(no,nb)]);
CB=CBT';% bar connectivity matrix
CR=CRT';% gravity center connectivity matrix of bars
         
% 2) for strings
% q00,...,q31 - sub connectivity matrices
e0=sparse(zeros(3,1));ee=speye(3);e1=ee(:,1);e2=ee(:,2);e3=ee(:,3);
q00=[-e1 e2 -e2 e3 -e3 e1];q01=sparse([1 -1 0 0 0 0]);q02=sparse([0 0 1 -1 0 0]);q03=sparse([0 0 0 0 1 -1]);
q11=kron(q01,e1);q12=kron(q01,e2);q13=kron(q01,e3);
q21=kron(q02,e1);q22=kron(q02,e2);q23=kron(q02,e3);
q31=kron(q03,e1);q32=kron(q03,e2);q33=kron(q03,e3);

% connectivity matrix for top and bottom strings 
if q==1 % for complexity q=1
    CStT=[q00;sparse(zeros(3,6));q12+q23+q31;sparse(zeros(3,6))];
    CSbT=[sparse(zeros(3,6));q00;sparse(zeros(3,6));q11+q22+q33];
    CSt=CStT';CSb=CSbT';
else % for complexity q>1
Iq=speye(q);Iq1=speye(q-1);
E1=sparse(1:q-1,2:q,ones(1,q-1),q,q);E2=sparse(2:q,1:q-1,ones(1,q-1),q,q);
E3=sparse(1:q-2,2:q-1,ones(1,q-2),q-1,q-1);E4=sparse(2:q-1,1:q-2,ones(1,q-2),q-1,q-1);
E5=sparse(1:q-2,1:q-2,ones(1,q-2),q-1,q-1);E6=sparse(3:q,1:q-2,ones(1,q-2),q*(q-1),q-1);
E7=sparse(2*q-2:-1:q+1,1:q-2,ones(1,q-2),2*q-1,q-2);E8=sparse(1:q-1,2:q,ones(1,q-1),q-1,q);


% connectivity matrix for top strings
psit00=q00;psit01=[sparse(zeros(3,6));q13;sparse(zeros(3*q*(q-1)-6,6))];psit02=[sparse(zeros(3,6));q23;sparse(zeros(3*q*(q-1)-6,6))];psit03=[sparse(zeros(3,6));q33;sparse(zeros(3*q*(q-1)-6,6))];

kt0=kron(Iq,q00)+kron(E1,q32);kt1=kron(E2,q13);kt2=kron(Iq,q21);
kt3=kron(E8,[q23;sparse(zeros(3*(q-1),6))]);kt4=[q33 sparse(zeros(3,6*(q-1)))];kt5=[sparse(zeros(1,6*(q-1))) q01];
kt7=[q03 sparse(zeros(1,6*(q-1)))];kt6=[kron(Iq,q01);sparse(zeros(q-2,6*q));kt7];
psit10=[q21 sparse(zeros(3,6*q*(q-1)-6))];psit11=kron(Iq1,kt0)+kron(E3,kt2)+kron(E4,kt1);
psit12=[kt3 sparse(zeros(3*q*(q-1),6*q*(q-2)))];psit13=kron(E6,kt4);psit14=[kron(E7,kt5) kt6];

psit20=[q22 sparse(zeros(3,6*q*(q-1)-6))];psit21=psit13;psit22=psit11;psit23=psit12;psit24=psit14;

psit30=[q23 sparse(zeros(3,6*q*(q-1)-6))];psit31=psit12;psit32=psit13;psit33=psit11;psit34=psit14;

PSIT0=[psit00;psit01;psit02;psit03;sparse(zeros(3*(3*q*(q-1)+1)+6*(2*q-1),6))];
PSIT1=[psit10;psit11;psit12;psit13;sparse(zeros(3*(3*q*(q-1)+1)+0*(2*q-1),6*q*(q-1)));psit14;sparse(zeros(5*(2*q-1),6*q*(q-1)))];
PSIT2=[psit20;psit21;psit22;psit23;sparse(zeros(3*(3*q*(q-1)+1)+1*(2*q-1),6*q*(q-1)));psit24;sparse(zeros(4*(2*q-1),6*q*(q-1)))];
PSIT3=[psit30;psit31;psit32;psit33;sparse(zeros(3*(3*q*(q-1)+1)+2*(2*q-1),6*q*(q-1)));psit34;sparse(zeros(3*(2*q-1),6*q*(q-1)))];
CStT=[PSIT0 PSIT1 PSIT2 PSIT3];% CStT - transpose of top string connectivity matrix 
CSt=CStT';% CSt - top string connectivity matrix

% connectivity matrix for bottom strings
psib00=q00;psib01=[sparse(zeros(3,6));q21;sparse(zeros(3*q*(q-1)-6,6))];psib02=[sparse(zeros(3,6));q31;sparse(zeros(3*q*(q-1)-6,6))];psib03=[sparse(zeros(3,6));q11;sparse(zeros(3*q*(q-1)-6,6))];

kb0=kron(Iq,q00)+kron(E1,q13);kb1=kron(E2,q21);kb2=kron(Iq,q32);
kb3=kron(E8,[q31;sparse(zeros(3*(q-1),6))]);kb4=[q11 sparse(zeros(3,6*(q-1)))];kb5=[sparse(zeros(1,6*(q-1))) q02];
kb7=[q01 sparse(zeros(1,6*(q-1)))];kb6=[kron(Iq,q02);sparse(zeros(q-2,6*q));kb7];
psib10=[q32 sparse(zeros(3,6*q*(q-1)-6))];psib11=kron(Iq1,kb0)+kron(E3,kb2)+kron(E4,kb1);
psib12=[kb3 sparse(zeros(3*q*(q-1),6*q*(q-2)))];psib13=kron(E6,kb4);psib14=[kron(E7,kb5) kb6];

psib20=[q33 sparse(zeros(3,6*q*(q-1)-6))];psib21=psib13;psib22=psib11;psib23=psib12;psib24=psib14;

psib30=[q31 sparse(zeros(3,6*q*(q-1)-6))];psib31=psib12;psib32=psib13;psib33=psib11;psib34=psib14;

PSIB0=[sparse(zeros(3*(3*q*(q-1)+1),6));psib00;psib01;psib02;psib03;sparse(zeros(6*(2*q-1),6))];
PSIB1=[sparse(zeros(3*(3*q*(q-1)+1),6*q*(q-1)));psib10;psib11;psib12;psib13;sparse(zeros(3*(2*q-1),6*q*(q-1)));psib14;sparse(zeros(2*(2*q-1),6*q*(q-1)))];
PSIB2=[sparse(zeros(3*(3*q*(q-1)+1),6*q*(q-1)));psib20;psib21;psib22;psib23;sparse(zeros(4*(2*q-1),6*q*(q-1)));psib24;sparse(zeros(1*(2*q-1),6*q*(q-1)))];
PSIB3=[sparse(zeros(3*(3*q*(q-1)+1),6*q*(q-1)));psib30;psib31;psib32;psib33;sparse(zeros(5*(2*q-1),6*q*(q-1)));psib34;sparse(zeros(0*(2*q-1),6*q*(q-1)))];
CSbT=[PSIB0 PSIB1 PSIB2 PSIB3];% CSbT - transpose of bottom string connectivity matrix 
CSb=CSbT';% bottom string connectivity matrix
end

% connectivity matrix for vertical strings
CSvT=sparse([-eye(nb);kron(eye(nu),[e2 e3 e1]);zeros(no,nb)]);% transpose of vertical string connectivity matrix
CSv=CSvT';% vertical string connectivity matrix

CST=[CStT CSbT CSvT];% transpose of string connectivity matrix
CS=CST';% string connectivity matrix
end

%*******************************
% 2. node matrix
%********************************
function [Node_I]=nodematrix(q,Height,r0)
% 1)vector dc1 and dc2
dc1=3*sqrt(2-sqrt(3))*r0*[-1/2;-sqrt(3)/2;0];dc2=3*sqrt(2-sqrt(3))*r0*[1;0;0];


% 2)calculate the center point of each unit Node_C(part0 and 1)
nc=1+q*(q-1);%number of center point(just consider part0 and 1)
Node_C=zeros(3,nc);
Node_C(:,1)=zeros(3,1);
if q>1
    for i=0:q-2
        for k=0:q-1
            n0=1+i*q+k+1;
            Node_C(:,n0)=(i+1)*dc1+k*dc2;
        end
    end
end


% 3)calculate the nodes of each unit
% the basic unit
NbasicT=zeros(3,3);NbasicB=zeros(3,3);
NbasicT(:,1)=r0*[cos(-pi/12-pi/2);sin(-pi/12-pi/2);Height/r0];NbasicT(:,2)=r0*[cos(7*pi/12-pi/2);sin(7*pi/12-pi/2);Height/r0];NbasicT(:,3)=r0*[cos(15*pi/12-pi/2);sin(15*pi/12-pi/2);Height/r0];
NbasicB(:,1)=r0*[cos(13*pi/12-pi/2);sin(13*pi/12-pi/2);0];NbasicB(:,2)=r0*[cos(21*pi/12-pi/2);sin(21*pi/12-pi/2);0];NbasicB(:,3)=r0*[cos(5*pi/12-pi/2);sin(5*pi/12-pi/2);0];

%bar nodes
%for part0 and 1
nu=1+3*q*(q-1);%number of unit
Node_Ib=zeros(3,6*nu);%part0
Node_Ib(:,1:3*nc)=kron(Node_C,ones(1,3))+kron(ones(1,nc),NbasicT);% part1 top nodes
Node_Ib(:,3*nu+1:3*nu+3*nc)=kron(Node_C,ones(1,3))+kron(ones(1,nc),NbasicB);% part1 bottom nodes

% for (2,i,j) and (3,i,j)
c12=cos(2*pi/3);s12=sin(2*pi/3);
C312=[c12 s12 0;-s12 c12 0;0 0 1];
c13=cos(4*pi/3);s13=sin(4*pi/3);
C313=[c13 s13 0;-s13 c13 0;0 0 1];
n11=4:3*(1+(q-1)*q);n12=3*nu+4:3*nu+3*(1+(q-1)*q);
n21=3*(1+(q-1)*q)+1:3*(1+2*(q-1)*q);n22=3*nu+3*(1+(q-1)*q)+1:3*nu+3*(1+2*(q-1)*q);
n31=3*(1+2*(q-1)*q)+1:3*(1+3*(q-1)*q);n32=3*nu+3*(1+2*(q-1)*q)+1:3*nu+3*(1+3*(q-1)*q);
Node_Ib(:,[n21 n22])=C312'*Node_Ib(:,[n11 n12]);% part2 bar nodes 
Node_Ib(:,[n31 n32])=C313'*Node_Ib(:,[n11 n12]);% part3 bar nodes

%string nodes
Node_Io=zeros(3,(6*q-3)*2);
%part 1
for i=1:q
    Node_Io(:,i)=Node_Ib(:,3)+q*dc1+i*dc2;
    Node_Io(:,i+6*q-3)=Node_Ib(:,1+3*nu)+q*dc1+i*dc2;
end
for i=1:q-2
    Node_Io(:,i+q)=Node_Ib(:,3)+(q-i)*dc1+q*dc2;
    Node_Io(:,i+q+6*q-3)=Node_Ib(:,1+3*nu)+(q-i)*dc1+q*dc2;
end
Node_Io(:,2*q-1)=Node_Ib(:,2)+(q-1)*dc1-1*dc2;
Node_Io(:,2*q-1+6*q-3)=Node_Ib(:,3+3*nu)+(q-1)*dc1-1*dc2;
no11=1:2*q-1;no12=3*(2*q-1)+1:4*(2*q-1);
no21=2*q-1+1:2*(2*q-1);no22=4*(2*q-1)+1:5*(2*q-1);
no31=2*(2*q-1)+1:3*(2*q-1);no32=5*(2*q-1)+1:6*(2*q-1);
Node_Io(:,[no21 no22])=C312'*Node_Io(:,[no11 no12]);%part2 string nodes
Node_Io(:,[no31 no32])=C313'*Node_Io(:,[no11 no12]);%part3 string nodes

Node_I=[Node_Ib Node_Io];% node matrix

end

