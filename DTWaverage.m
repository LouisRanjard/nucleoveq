function [ dist, mat4, refpos ] = DTWaverage( mat1, mat2, cow, q, ce, fe )
% perform a dynamic time warping in order to measure a distance between matrices
% according to the algorithm of Oommen1995
% Kruskal&Liberman in Chapt4 Sankoff&Kruskal1999 (time warps...)
% combining insertion/deletion and compression/expansion
% if s_ave==1 then it computes the average between the 2 vector sequences
% q is the weight for the computation of the average (must be between 0 and 1)
%        it affects mat1 and mat2 has the weight 1-q
% cow is the weight of each coefficient for the distance computation, must be the same length than size(mat1,1)
% ce indicates whether to use compression/expansion (1) or not (0)
% fe indicates whether to use free-end alignment (1) or not (0)
% an average sequence can be computed:
    % bak=1 expansion
    % bak=2 insertion
    % bak=3 substitution
    % bak=4 deletion
    % bak=5 compression
%
% dist: alignment distance
% mat4: average vector sequence 
% refpos: position where alignment FINISHES in mat1
%

% if a matrix is empty return -1
if size(mat1,2)==0 || size(mat2,2)==0
    dist=-1;
    mat4=[];
    % error('DTWaverage(): at least of the input matrices is empty');
    return;
end

switched = 0;
if ( size(mat1,1)~=size(mat2,1) && (nargin<6 || fe==0) ) % fe: not relevant when aligning to reference (persistence row)
    if size(mat1,2)==size(mat2,2)
        mat1=mat1';
        mat2=mat2';
        switched = 1;
    else
        error('DTWaverage(): matrix sizes are not compatible');
    end
end

% initialise the coefficients weight by one if nothing specify or if cow==1
if (nargin<3 || numel(cow)==0)
    % all the same weight
    cow=ones(1,size(mat1,1));
    % for 12PLPlogE tieke
    %cow=[0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4];
elseif cow==1
    cow=ones(1,size(mat1,1));
end

% if no sequence weight specify do not compute average sequence
if (nargin<6)
    fe = 0 ;
    if (nargin<5)
        ce = 0 ;
        if (nargin<4)
            q = 0 ;
        end
    end
end

if (fe && size(mat1,1)~=5 && size(mat1,1)~=4)
    error('DTWaverage(): Matrix dimension issue');
end

%%%%%%%%%%%%%% CALL THE C FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% need to be compiled first: mex GCC=/usr/bin/gcc-4.7 COPTIMFLAGS=-O3 [path]/DTWave.c
[dist, mat4, refpos] = DTWave(mat1,mat2,cow,q,ce,fe);

% return the matrix in the same direction as input
if switched==1, mat4 = mat4'; end

% debugging, controls
%fprintf(1,'%d\n',refpos);
% if sum(sum(round(mat1(1:4,:).*10000)/10000,1)~=1) || sum(sum(mat2(1:4,:),1)~=1) 
%     fprintf(1,'error\n');
% end
