% 358 Project 1
% Matthew Ng
% Most of the code was given by Prof. Murali Subbarao
% 
% Summary:
% Displays a cube that roates and falls.
% One face contains and image. 




%define 8 points of the cube in world coordinate
length = 10;
V1= [ 0 0 0 ];
V2= [ 0 length 0 ];
V3= [ length length 0 ];
V4= [ length 0 0 ];
V5= [ length 0 length ];
V6= [ 0 length length ];
V7= [ 0 0 length ];
V8= [ length length length ];

% Find the unit vector u81 corresponding to the axis of rotation which 
%is along (V8-V1).From u81, compute the 3x3 matrix N in Eq. 2.32 used for
%computing the rotation matrix R in eq. 2.34
%?????????????????
v81 = V8-V1;
mag81 = norm(v81);

u81 = v81/mag81;

Nx = [0,- u81(3), u81(2); u81(3), 0, - u81(1) ; - u81(2) , u81(1) ,0];

T0 = [ -20 -25 500 ]; % origin of object coordinate system in mm
%T0 = [ -30 -20 500 ]; % origin of object coordinate system

%set given values
f=40; %focal length in mm
%f=input(&#39;f=&#39;);
%Initialize the 3x3 camera matrix K given the focal length
K= [f,0,0 ; 0,f,0 ; 0,0,1]; %?????
velocity=[2 9 7 ]; % translational velocity
theta0=0;
w0=20;% angular velocity in deg/sec
p=0.01;%pixel size(mm)
Rows = 600; %image size
Cols = 600; % image size
A = zeros(Rows, Cols); %output image
%B = rgb2gray(imrotate(imresize(imread('3Dbg.jpg'), [Rows,Cols]), -90));
B = imread('stadium1.jpg');
%comment the line above and the line 125 to get black background
r0= round(Rows/2);
c0 = round(Cols/2);

% You are given a rectangle/square in 3D space specified by its
% corners at 3D position vectors V1, V2, V3, V4.
% You are also given a rectangular/square graylevel image
% tmap of size r x c.
% This image is to be 'painted' on the 3D rectangle/square, and
% for each pixel at position (i,j),
% the corresponding 3D coordinates
% X(i,j), Y(i,j), and Z(i,j), should be computed,
% and that 3D point is
% associated with the brightness given by tmap(i,j).
%
% Find the unit vectors corresponding to u21=(V2-V1)/|(V2-V1)|

% and u41= (V4-V1)/|(v4-V1), and compute X(i,j), Y(i,j), and
%Z(i,j).
% Compute the unit vector u21 along (V2-V1) and
% Compute the unit vector u41 along (V4-V1) and
h= norm(V2 -V1); % height = distance from v2 to v1
w= norm(V4 - V1); % width = distance from v4 to v1
u21= (V2-V1)/h;
u41= (V4-V1)/w;

% For each pixel of texture map, compute its (X,Y,Z) values
%
tmap = imread('einstein50x50v.jpg'); % texture map image
[ r c ] = size(tmap);
X=zeros(r,c);
Y=zeros(r,c);
Z=zeros(r,c);
for i = 1 : r
    for j = 1 : c 
        p1 = V1 + (i-1)* u21* (h/r) + (j-1)* u41*(w/c);
        X(i,j)=p1(1);
        Y(i,j)= p1(2);%filled this in
        Z(i,j)= p1(3);%filled this in
    end
end
acc = [ 0.0 -0.80 0 ]; %acceleration
for t=0:0.2:24 % Generate a sequence of images
                % as a function of time
    theta=theta0+w0*t;
    T=T0+ velocity*t + 0.5 * acc * t * t;
    % Compute the Rotation matrix from N and theta
    R= eye(3) + sind(theta)*Nx +((1-cosd(theta))*Nx*Nx); %????????????????
    %find the image position of vertices
    v=Map2Da(K,R,T,V1);
    v1 = MapIndex(v,c0,r0,p);
    
    v=Map2Da(K,R,T,V2);
    v2 = MapIndex(v,c0,r0,p);
    
    v=Map2Da(K,R,T,V3);
    v3 = MapIndex(v,c0,r0,p);
    
    v=Map2Da(K,R,T,V4);
    v4 = MapIndex(v,c0,r0,p);
    
    v=Map2Da(K,R,T,V5);
    v5 = MapIndex(v,c0,r0,p);
    
    v=Map2Da(K,R,T,V6);
    v6 = MapIndex(v,c0,r0,p);
    
    v=Map2Da(K,R,T,V7);
    v7 = MapIndex(v,c0,r0,p);
    
    v=Map2Da(K,R,T,V8);
    v8 = MapIndex(v,c0,r0,p);

    %?????????????????????????????
    % Draw edges of the cube
    A = zeros(Rows, Cols);
    A = B;%comment this line out to get black background
    A = Line(A, v1,v2);
    %????????????????????????????
    A = Line(A,v2,v3);
    A = Line(A,v3,v4);
    A = Line(A,v4,v1);
    A = Line(A,v7,v6);
    A = Line(A,v6,v8);
    A = Line(A,v8,v5);
    A = Line(A,v5,v7);
    A = Line(A,v1,v7);
    A = Line(A,v2,v6);
    A = Line(A,v3,v8);
    A = Line(A,v4,v5);
    
    

    % Add texture map to one face.
    for i = 1 : r
        for j = 1 : c
            p1 = [ X(i,j) Y(i,j) Z(i,j) ];             
% %Find the row and column indices (ir,jr) in integers that give
% %the image position of point p1 in A.
% % Use the same method as for the corners of the cube above.
% 
% %?????????????????????????
            v = Map2Da(K,R,T,p1);
            mapped = MapIndex(v,c0,r0,p);
            ir = round(mapped(1));
            jr = round(mapped(2));
%             disp('ir = ')
%             disp(ir)
%             disp('jr = ')
%             disp(jr)
            %breakpoint
            if((ir>0) && (jr>0) && (ir<=Rows) && (jr<=Cols))               
                A(ir,jr)=tmap(i,j);
%                 disp("A =")
%                 disp(A(ir,jr))
%                 disp("tmap =")
%                 disp(tmap(i,j))
            end
        end
% In a general case, you may need to
% fill up gaps in A(ir,jr)
% through interpolation. But, in this project,
% you can skip interpolation. The output will not
% look nice due to the gaps.
    end

    A=mat2gray(A);
    imshow(A);
    %pause
% pause if you want to display frame by frame
% and press return to display the next frame
end

%function for rotation and translation
function [ v ] =Map2Da( K,R,T,Vi)
    P=K*[R T']*[Vi 1]';
    %disp("p = ")
    %disp(P)
    w1=P(3,1);
    v(1)=P(1,1)/w1;
    %disp(v)
    v(2)= P(2,1)/w1;%?????????????;
end

%function for mapping image coordinates in mm to
% row and column index of the image, with pixel size p mm and
% image center at (r0,c0)
function [ v ] =MapIndex( u,c0,r0,p )
    v(1)= round(r0-u(2)/p);
    v(2)=round(c0+u(1)/p);
end
%function for line
% Draw line from v1 to v2 in image A
function [A] = Line(A, v1, v2)
    d=sqrt((v1-v2)*(v1-v2)');
    ui=(v2(1)-v1(1))/d;
    uj=(v2(2)-v1(2))/d;
    i=v1(1);
    j=v1(2);
    [ rows cols ] = size(A);
    for K=0:round(d)
        i=i+ui;
        j=j+uj;
        ir=round(i);
        jr=round(j);
        if((ir>0)&&(jr>0) && (ir<rows) && (jr<cols))
            A(ir,jr)=255;
        end
    end
end