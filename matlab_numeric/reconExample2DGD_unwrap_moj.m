% very basic and crude non-Cartesian recon using griddata()
%
% needs mapVBVD in the path

%% Load the latest file from a dir
%path='/ad/O/Exchange/Zaitsev/AMR_LFOV';
path = 'C:\Users\shafiekh\Documents\MATLAB\MojiMAT\tutorials-main\11_from_GRE_to_EPI\data';
%path = '/home/wehkamp/netzlaufwerke/O/Exchange/Zaitsev';

%path='../IceNIH_RawSend/'; % directory to be scanned for data files
%path='~/Dropbox/shared/data/siemens/';
%path='~/Downloads/2021-07-12-090810/';
pattern='*.seq';
nF = 13;
if path(end)~=filesep, path(end+1)=filesep; end
D=dir([path pattern]);
%[~,I]=sort([D(:).datenum]);
[~,I]=sort(string({D(:).name}));
seq_file_path = [path D(I(nF)).name]; % use end-1 to reconstruct the second-last data set, etc.

%% alternatively just provide the path to the .seq file
%seq_file_path='../interpreters/siemens/data_example/gre_example.seq'
%seq_file_path='/data/20211025-AMR/data/2021-10-26-093137.seq';

%% Load sequence from the file with the same name as the raw data
fprintf(['loading `' seq_file_path '´ ...\n']);
seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');

%% keep basic filename without the extension
[p,n,e] = fileparts(seq_file_path);
basic_file_path=fullfile(p,n);

%% load the raw data file
data_file_path= [basic_file_path '.mat']; % try to load a matlab file with raw data...
try
    fprintf(['loading `' data_file_path '´ ...\n']);
    data_unsorted = load(data_file_path);
    if isstruct(data_unsorted)
        fn=fieldnames(data_unsorted);
        assert(length(fn)==1); % we only expect a single variable
        data_unsorted=double(data_unsorted.(fn{1}));
    end
catch
    data_file_path= [basic_file_path '.dat']; % now try to load a raw data file...
    fprintf(['falling back to `' data_file_path '´ ...\n']);
    twix_obj = mapVBVD(data_file_path);
    if iscell(twix_obj)
        data_unsorted = double(twix_obj{end}.image.unsorted());
    else
        data_unsorted = double(twix_obj.image.unsorted());
    end
    seqHash_twix=twix_obj.hdr.Dicom.tSequenceVariant;
    if length(seqHash_twix)==32
        fprintf(['raw data contain pulseq-file signature ' seqHash_twix '\n']);
    end
    %clear twix_obj
end

%% calculate k-space trajectory
%[3.5 4 0] % for AMR 2D-UTE
traj_recon_delay=0*1e-6;% [0.527 -1.367 0]; % adjust this parameter to potentially improve resolution & geometric accuracy. It can be calibrated by inverting the spiral revolution dimension and making two images match. for our Prisma and a particular trajectory we found 1.75e-6
grad_offsets=[0 0 0];

seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');
%[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace('trajectory_delay', traj_recon_delay);
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay, 'gradient_offset', grad_offsets);
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); 
title('2D k-space trajectory'); drawnow;

%% Define FOV and resolution and simple off-resonance frequency correction 

fov=650e-3; Nx=256; Ny=Nx; % define parameters explicitly
Ns=1; % this is a number of slices (or contrasts) which needs to be specified manually for now
%Na=1;

def_fov=seq.getDefinition('FOV'); % try to read from definitions 
if numel(def_fov)
    fov=max(def_fov);
end
%fov=fov*1.5;
deltak=1/fov;

% or estimate from the k-space trajectory
%k_max=max(vecnorm(ktraj_adc));
%Nx=round(k_max/deltak*2);
%Ny=Nx; 

os=2; % oversampling factor (we oversample both in image and k-space)
offresonance=500; % global off-resonance in Hz

%%
rawdata = permute(data_unsorted, [1,3,2]);
adc_len=size(rawdata,1);
readouts=size(rawdata,2);
channels=size(rawdata,3);

Na=numel(rawdata(:,:,1))/Ns/numel(t_adc); % averages/acquisitions
if Na>1
    nTRs=size(rawdata,2)/Na;
    rawdata = sum(permute(reshape(rawdata, [size(rawdata,1),size(rawdata,2)/Na,Na,size(rawdata,3)]),[1,2,4,3]),4);    
    readouts=readouts/Na;
end

if strcmp('ute_rs',seq.getDefinition('Name'))
    % average every 2nd spoke because of the half-rf excitation
    ktraj_adc=reshape(ktraj_adc,[3,adc_len,readouts]);
    ktraj_adc=ktraj_adc(:,:,1:2:end-1);
    t_adc=reshape(t_adc,[1,adc_len,readouts]);
    t_adc=t_adc(:,:,1:2:end-1);
    %rawdata=rawdata(:,1:2:end-1,:);
    rawdata=rawdata(:,1:2:end-1,:)+rawdata(:,2:2:end,:);
    readouts=readouts/2;
    ktraj_adc=reshape(ktraj_adc,[3,adc_len*readouts]);
    t_adc=reshape(t_adc,[1,adc_len*readouts]);
end

rawdata = reshape(rawdata, [size(rawdata,1)*size(rawdata,2)/Ns,Ns,size(rawdata,3)]);
ktraj_adc=ktraj_adc(:,1:end/Ns);
t_adc=t_adc(1:end/Ns);

t_adcEx=t_adc-t_excitation(end);
for i=1:length(t_excitation)-1
    ii=find(t_adc>t_excitation(i)&t_adc<t_excitation(i+1));
    t_adcEx(ii)=t_adc(ii)-t_excitation(i);
end

for s=1:Ns
    for c=1:channels
        rawdata(:,s,c) = rawdata(:,s,c) .* exp(-1i*2*pi*t_adcEx'*offresonance);
    end
end

%% here we expect Nx, Ny, deltak to be set already
% and rawdata ktraj_adc loaded (and having the same dimensions)

kxm=round(os*os*Nx/2);
kym=round(os*os*Ny/2);

[kyy,kxx] = meshgrid(-kxm:(kxm-1), -kym:(kym-1));
kyy=-kyy*deltak/os; % we swap the order and invert one sign to account for Matlab's strange column/line convention
kxx=-kxx*deltak/os; % this is needed to match Localizer images on TRIO

kgd=zeros([size(kxx) Ns channels]);
d1=1; d2=size(rawdata,1);
%d1=320*64+1; d2=320*192;
%d1=1; d2=320*128;
%d1=320*128+1; d2=320*256;
tic;
for s=1:Ns
    for c=1:channels        
        kgd(:,:,s,c)=griddata(ktraj_adc(1,d1:d2),ktraj_adc(2,d1:d2),rawdata(d1:d2,s,c),kxx,kyy,'cubic');
        %F=scatteredInterpolant(ktraj_adc(1,d1:d2)',ktraj_adc(2,d1:d2)',squeeze(rawdata(d1:d2,s,c)),'natural','none');     
        %kgd(:,:,s,c)=F(kxx,kyy);
    end
end
kgd(isnan(kgd))=0;
toc

figure;imagesc(log(abs(kgd(:,:,1,1))));axis('square');
title('2D k-space data, ch1 as log(abs())');

igd=ifftshift(ifft2(ifftshift(kgd)));

Nxo=round(Nx*os);
Nyo=round(Ny*os);
Nxs=round((size(igd,1)-Nxo)/2);
Nys=round((size(igd,2)-Nyo)/2);
images = igd((Nxs+1):(Nxs+Nxo),(Nys+1):(Nys+Nyo),:,:);

%% Image display with optional sum of squares combination
figure;
if channels>1
    sos=abs(sum(images.*conj(images),ndims(images))).^(1/2);
    sos=sos./max(sos(:));    
    imab(sos); title('reconstructed image(s), sum-of-squares');
    %imwrite(sos, ['img_combined.png']
else
    imab(abs(images)); title('reconstructed image(s)');
end
colormap('gray');
saveas(gcf,[basic_file_path '_image_2d_gridding'],'png');

%% create voxel geometry struct 

fprintf('ATTENTION, untested code ahead!\n');

% image size
vg1.len1=size(sos,1);
vg1.len2=size(sos,2);
vg1.len3=1; 

% image resolution
vg1.vsize1 = 1000/vg1.len1/deltak;%twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV * ro_os / vg.len1; % alternative: 1/k1_max
vg1.vsize2 = 1000/vg1.len2/deltak;%twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV / vg.len2; % alternative: 1/k2_max
vg1.vsize3 = 1;%twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness / 16; %vg.depth; % alternative: 1/k3_max; 

% % volume orientation 
qv=twix_obj.image.slicePos(4:7,1);
tv=twix_obj.image.slicePos(1:3,1);
q=quaternion(qv');
rm=quat2rotm(q)';

% signs copied from freesurfer, probably related to the patient orientation
vg1.v2_ras = [-1 -1 +1].*rm(1,:);
vg1.v1_ras = [-1 -1 +1].*rm(2,:);
vg1.v3_ras = [-1 -1 +1].*rm(3,:);

% volume offset
vg1.c_ras = [-1 -1 +1].*tv';

% vg_i_to_r


%% calculate deformation field
gradUnwarp = GradUnwarpV();

gradUnwarp.read_siemens_coeff('~/range_software/pulseq/matlab/idea/grad/coeff_AS82.grad');
%gradUnwarp.read_siemens_coeff('C:\Users\shafiekh\Desktop\gradient_nonlinearity_correction-main (1)\gradient_nonlinearity_correction-main\matlab\coeff_AS82.grad');

%gradUnwarp.read_siemens_coeff('coeff_GC98SQ.grad');

% x0 = struct2array(load("x0_new5","-mat"));
% gradUnwarp.update_AB_from_sph_coeff(x0);

tic
gradUnwarp.create_transtable(vg1);
toc

%% unwarp!
%figure; imagesc(sos(end:-1:1,end:-1:1,ceil(end/2))); title('warped slice');
figure; imagesc(rot90(sos(end-5:-1:8,end-4:-1:9,ceil(end/2)),2)); title('Warped Slice'); colorbar

F=griddedInterpolant({0:(vg1.len1-1),0:(vg1.len2-1)},sos(end:-1:1,end:-1:1),'cubic','none');
V1u = gradUnwarp.jd.*F(gradUnwarp.gcam(:,:,:,1),gradUnwarp.gcam(:,:,:,2));

clear F;
V1u(V1u<0)=0;
%Opt_based_unwarp = V1u;
Siemens_based_unwarp = V1u;

%figure; imagesc(V1u(:,:,ceil(end/2))); title('unwarped slice');
figure; imagesc(V1u(end-8:-1:5,end-8:-1:5,ceil(end/2))); title('Unwarped Slice using Optimization-based Coefficients'); colorbar
%figure; imagesc(V1u(end-8:-1:5,end-8:-1:5,ceil(end/2))); title('Unwarped Slice using Siemens Coefficients'); colorbar
caxis([0,max(sos(:))]);

figure();
subplot(131);imagesc(rot90(sos(end-5:-1:8,end-4:-1:9,ceil(end/2)),2)); title('Warped Slice'); colorbar
subplot(132);imagesc(Siemens_based_unwarp(end-8:-1:5,end-8:-1:5,ceil(end/2))); title('Unwarped Slice using Siemens Coefficients'); colorbar
subplot(133);imagesc(Opt_based_unwarp(end-8:-1:5,end-8:-1:5,ceil(end/2))); title('Unwarped Slice using Optimization-based Coefficients'); colorbar
