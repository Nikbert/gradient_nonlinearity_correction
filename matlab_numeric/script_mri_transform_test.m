%% volume geometry for the entire FOV
R0=250; %Radius in mm
dx=10;  %Diskretisierung
ns=length(-R0:dx:R0);

% image size
vg.len1=ns;
vg.len2=ns;
vg.len3=ns;

% image resolution
vg.vsize1 = dx;
vg.vsize2 = dx;
vg.vsize3 = dx; 

% volume orientation 
vg.v1_ras = [1 0 0];
vg.v2_ras = [0 1 0];
vg.v3_ras = [0 0 1];

% volume offset
vg.c_ras = [0 0 0];

%%
gradUnwarp = GradUnwarpV();

%gradUnwarp.read_siemens_coeff('coeff_GC98SQ.grad');
% gradUnwarp.read_siemens_coeff('~/range_software/pulseq/matlab/idea/grad/coeff_AS82.grad');
gradUnwarp.read_siemens_coeff('coeff_AS82.grad');

tic
gradUnwarp.create_transtable(vg);
toc

%% visualize distortion maps
% this is relatively computationally expensive visualization, but nice...
xx=zeros(3,vg.len1,vg.len2,vg.len3);
[xx(1,:,:,:),xx(2,:,:,:),xx(3,:,:,:)]=ndgrid(0:(vg.len1-1),0:(vg.len2-1),0:(vg.len3-1));
dd=permute(gradUnwarp.gcam,[4,1,2,3])-xx;
figure; imagesc(squeeze(vecnorm(dd(:,:,:,ceil(end/2))))); title('distortion distance, central transverse slice'); colorbar;
axis('square'); % enforce aspect ratio
figure; imagesc(squeeze(vecnorm(dd(:,:,ceil(end/2),:)))); title('distortion distance, central coronal? slice'); colorbar;
axis('square'); % enforce aspect ratio

% amplitude scaling visualisation 
figure; imagesc(gradUnwarp.jd(:,:,ceil(end/2))); title('Jacobian determinant (amplitude scaling)'); colorbar;
axis('square'); % enforce aspect ratio
figure; imagesc(squeeze(gradUnwarp.jd(:,ceil(end/2),:))); title('Jacobian determinant (amplitude scaling), coronal?'); colorbar;
axis('square'); % enforce aspect ratio

%% now create the distortion map directly in device coordinates
[xx(1,:,:,:),xx(2,:,:,:),xx(3,:,:,:)]=ndgrid(-R0:dx:R0,-R0:dx:R0,-R0:dx:R0);
iSphere=find(vecnorm(xx(:,:))<=R0);

fprintf('creating distortion map\n');
tic;
xxd=reshape(gradUnwarp.lab_to_mri(xx(:,:)),size(xx));
toc

dd=xxd-xx;
figure; imagesc(squeeze(vecnorm(dd(:,:,:,ceil(end/2))))); title('distortion distance (mm), central transverse slice'); colorbar;
axis('square'); % enforce aspect ratio
figure; imagesc(squeeze(vecnorm(dd(:,:,ceil(end/2),:)))); title('distortion distance (mm), central coronal? slice'); colorbar;
axis('square'); % enforce aspect ratio

%% create a basis set for fitting

points=xx(:,:);%xx(:,iSphere);

nmax=25;

fprintf('creating sh basis\n');
tic;
[basis, basis_normalization, ns, ms, ab]=create_sph_basis(points,R0,nmax);
toc
nbasis=size(basis,1);

%% convert Siemens coeficients
coefS3=zeros(nbasis,3);
for i=1:size(gradUnwarp.nm,2)
    qA=find(ns==gradUnwarp.nm(1,i)&ms==gradUnwarp.nm(2,i)&ab==1);
    qB=find(ns==gradUnwarp.nm(1,i)&ms==gradUnwarp.nm(2,i)&ab==2);
    for j=1:3
        coefS3(qA,j)=gradUnwarp.A3(j,i);
        coefS3(qB,j)=gradUnwarp.B3(j,i);
    end
end

coefS3(nmax+2,1)=-1; %x
coefS3(nmax+3,2)=-1; %y
coefS3(2,3)=1; %z

%% fit!
Bz=squeeze(xxd(1,:,:,:))/R0;
tic;

%coef=basis(:,iSphere)'\Bz(iSphere)';
basis_pinv=pinv(basis(:,iSphere)');
coef=basis_pinv*Bz(iSphere)';
toc

%BzRes=reshape(basis'*coef,size(Bz));
BzRes=zeros(size(Bz))+NaN;
BzRes(iSphere)=basis(:,iSphere)'*coef;

figure;imagesc(BzRes(:,:,ceil(end/2)));colorbar;
figure;imagesc(squeeze(BzRes(:,ceil(end/2),:)));colorbar;
figure;imagesc(squeeze(BzRes(:,ceil(end/2),:)-Bz(:,ceil(end/2),:)));colorbar;
figure;imagesc(squeeze(BzRes(ceil(end/2),:,:)));colorbar;
fprintf('maximum deviation in the sphere: %g\n',max(max(max(abs(Bz-BzRes)))));

%% compare with Siemens

figure;plot(coef);
hold on;plot(coefS3(:,1).*basis_normalization,'o'); 

%% create correct reverse mapping
fprintf('creating accurate reverse distortion map\n');
tic;
xxr=reshape(gradUnwarp.mri_to_lab(xx(:,:)),size(xx));
toc

%save

%% visualize the fields
iSphereDist=find(vecnorm(xxr(:,:))<=R0);

maskSphere=zeros(size(squeeze(xx(1,:,:,:))));
maskSphere(iSphere)=1;
maskSphereDist=zeros(size(squeeze(xx(1,:,:,:))));
maskSphereDist(iSphereDist)=1;

figure; contour(squeeze(maskSphere(:,ceil(end/2),:)),[0.5 0.55]);
hold on; contour(squeeze(xxd(1,:,ceil(end/2),:)),-R0:5:R0);
hold on; contour(squeeze(xxd(3,:,ceil(end/2),:)),-R0:5:R0);
axis('square'); title('distortins and FOV in laboratory coordinates');

figure; contour(squeeze(maskSphereDist(:,ceil(end/2),:)),[0.5 0.55]);
hold on; contour(squeeze(xxr(1,:,ceil(end/2),:)),-R0:5:R0);
hold on; contour(squeeze(xxr(3,:,ceil(end/2),:)),-R0:5:R0);
axis('square'); title('distortins and FOV in distorted MRI coordinates');

%% fit accurate reverse mapping with harmonics

fprintf('fitting reverse distortion maps with the original basis\n');
BzRes=zeros(size(Bz))+NaN;

tic;
%basis_pinv=pinv(basis');
basis_pinv=pinv(basis(:,iSphereDist)');

for i=1:3
    Bz=squeeze(xxr(i,:,:,:))/R0;

    %coef_rev_acc{i}=basis'\Bz(iSphereDist)';
    coef_rev_acc{i}=basis_pinv*Bz(iSphereDist)';
    
    BzRes(iSphereDist)=basis(:,iSphereDist)'*coef_rev_acc{i};
    fprintf('channel %d: maximum deviation in the distorted sphere: %g\n',i,max(max(max(abs(Bz-BzRes)))));
end
toc

%%
figure; plot(coef_rev_acc{1});
hold on; plot(coef_rev_acc{2});
hold on; plot(coef_rev_acc{3});

%% create a basis set for reverse fitting

fprintf('creating reverse basis\n');
tic;
points=xxd(:,iSphereDist);
[basis_rev, basis_normalization_rev, ns_rev, ms_rev, ab_rev]=create_sph_basis(points,R0,nmax);
toc

%% fit linear functions in distorted coordinates
fprintf('fitting linears in reverse basis\n');
tic;
basis_rev_pinv=pinv(basis_rev');
for i=1:3
    Bz=squeeze(xx(i,:,:,:))/R0;

    %coef_rev{i}=basis_rev'\Bz(iSphereDist)';
    coef_rev{i}=basis_rev_pinv*Bz(iSphereDist)';
end
toc

figure; plot(coef_rev{1});
hold on; plot(coef_rev{2});
hold on; plot(coef_rev{3});

%% use reverse basis to predict reverse distortion map 
points=xx(:,iSphereDist);
xxr_syn=zeros(size(xx))+NaN;
for i=1:3
    xxr_syn(i,iSphereDist)=basis(:,iSphereDist)'*coef_rev{i};
    fprintf('channel %d: maximum deviation in the distorted sphere: %g\n',i,max(abs(xxr_syn(i,:)-xxr(i,:)/R0)));
end

%%
% %% unwarp!
% if size(V1,3)>1
%     F=griddedInterpolant({0:(vg.len1-1),0:(vg.len2-1),0:(vg.len3-1)},double(V1),'cubic','none');
%     V1u = gradUnwarp.jd.*F(gradUnwarp.gcam(:,:,:,1),gradUnwarp.gcam(:,:,:,2),gradUnwarp.gcam(:,:,:,3));
% else
%     F=griddedInterpolant({0:(vg.len1-1),0:(vg.len2-1)},double(V1),'cubic','none');
%     V1u = gradUnwarp.jd.*F(gradUnwarp.gcam(:,:,:,1),gradUnwarp.gcam(:,:,:,2));
% end
% 
% clear F;
% V1u(V1u<0)=0;
% figure; imagesc(V1u(:,:,ceil(end/2))); title('unwarped central slice');
% axis('square'); % enforce aspect ratio
% 
% %% load Siemens' correction result
% V2 = squeeze(dicomreadVolume('./dicom_3d_offcenter_dc'));
% %V2 = squeeze(dicomreadVolume('./dicom_3d_amr_dc'));
% %V2 = squeeze(dicomreadVolume('./dicom_3d_off_amr_dc'));
% figure; imagesc(V2(:,:,ceil(end/2))); title('siemens-unwarped central slice');
% axis('square'); % enforce aspect ratio
% 
% figure; imagesc(double(V2(:,:,ceil(end/2)))-V1u(:,:,ceil(end/2))); title('difference to siemens, central slice');colorbar;
% axis('square'); % enforce aspect ratio
% 
% %% load freesurfer's warpping table 
% [vol_orig, vol_dest, vol_ind0, spacing, exp_k] = mris_read_m3z('unwrap_3d_offcenter_table.m3z');
% %vol_orig_xd(:,:,:,1)=vol_orig(:,:,:,1)-xx;
% %figure;imagesc(vol_orig_xd(:,:,1));
% %vol_dest_xd(:,:,:,1)=vol_dest(:,:,:,1)-xx;
% %figure;imagesc(vol_dest_xd(:,:,end/2));
% figure;imagesc(vol_dest(:,:,end/2,1)-gradUnwarp.gcam(:,:,end/2,2)'); title('distortion map difference, X, central slice'); colorbar;
% axis('square'); % enforce aspect ratio
% 
% %% load freesurfer's correction result
% test=MRIread('unwrap_3d_offcenter.nii');
% figure;imagesc(test.vol(:,:,end/2)); title('freesurfer-unwarped central slice'); 
% axis('square'); % enforce aspect ratio
% 
% %% test MRI <-> LAB coordinate transforms
% gradUnwarp = GradUnwarpV();
% gradUnwarp.read_siemens_coeff('~/range_software/pulseq/matlab/idea/grad/coeff_AS82.grad');
% xd=gradUnwarp.lab_to_mri([[100;150;200] [-200;-150;-100]]);
% gradUnwarp.mri_to_lab(xd)

function [basis, basis_normalization, ns, ms, ab]=create_sph_basis(points,R0,nmax)

npoints=size(points,2);

RV = vecnorm(points);
cosThetaV = points(3,:)./RV;
cosThetaV(abs(RV)<eps)=1; % fix singularities
PhiV = atan2(points(2,:), points(1,:));

RV=RV/R0; % pre-normalize RV
        
nbasis=(nmax+2)*(nmax+1)-1-nmax;
basis=zeros(nbasis,npoints);
ns=zeros(1,nbasis);
ms=zeros(1,nbasis);
ab=zeros(1,nbasis);%a=1 b=2
c=1;
for m=0:nmax
    nn=m:nmax;
    Pa=gsl_sf_legendre_Plm_e_mult(nn, m, cosThetaV);
    %assert(all(max(abs(Pa'))));
    
    cosPhi = cos(m*PhiV); 
    sinPhi = sin(m*PhiV); 
            
    if m>1
        F = RV .^ (m-1);
    else
        F = 1;
    end
    s=0;
    a=1;
    for n=nn
        if n>0
            F = F .* RV; 
        end
        basis(c+s,:)=F.*Pa(a,:).*cosPhi;
        ns(c+s)=n;
        ms(c+s)=m;
        ab(c+s)=1;
        s=s+1;
        if m>0
            basis(c+s,:)=F.*Pa(a,:).*sinPhi;
            ns(c+s)=n;
            ms(c+s)=m;
            ab(c+s)=2;
            s=s+1;
        end
        a=a+1;
    end
    c=c+s;
end
assert(c-1==nbasis);

% normalize
basis_normalization=max(abs(basis),[],2);
assert(all(basis_normalization>0));
basis=basis./basis_normalization(:,ones(1,npoints));

end
