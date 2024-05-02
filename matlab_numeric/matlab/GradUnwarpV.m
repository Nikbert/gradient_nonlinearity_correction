classdef GradUnwarpV < handle
% gradient unwarping class, initially modelled from the corresponding 
% freesurfer's implementation, but converted to vector maths in MAtlab and 
% heavily optimized for performance by Maxim Zaitsev 

properties
    % input parameters
    R0_mm; % nominal radius
    nm;    % list of n/m pairs for which A and B coeficients are defined
    A3;    % list of the A-coefficients for the three gradient axes
    B3;    % list of the B-coefficients for the three gradient axes
    AB;    % compressed list of the A and B-coefficients for the three gradient axes
    coeff;
    coeffDim;
    AB_initialized;
    gcam; % distortion table
    jd;   % jacobian determinant (amplitude correction factor)
end
    
methods
    function obj = GradUnwarpV()
        obj.AB_initialized=false;
    end
        
    function DV = spharm_evaluate(obj, XV)
        % function that evaluates the valies of the spherial (solid)
        % harmonic expansion described by the initialized class for the
        % given vector of points in 3D space
        
        if ~obj.AB_initialized
            error("gradient file not loaded!");
        end
        
        % nm : list of used pairs of n amd m indexes
        % A3 and B3 : all Alphas and Betas
        
        % convert input points to spherical coordinates
        RV = vecnorm(XV);
        cosThetaV = XV(3,:)./RV;
        cosThetaV(abs(RV)<eps)=1; % fix simgularities
        PhiV = atan2(XV(2,:), XV(1,:));
        
        RV=RV/obj.R0_mm; % pre-normalize RV to save computations in the loop
        
        % initialize the output
        DV = zeros(size(XV));
        
        % the function gsl_sf_legendre_Plm_e_mult() returns harmonics for
        % all n's (degrees) for given m (order). Therefore the outer loop
        % eeds to be over unique m's, which saves about 50% of calculations
        mm=unique(obj.nm(2,:));
        for m=mm
            nn=unique(obj.nm(1,obj.nm(2,:)==m)); % now extract uniqe n's for the given m
            
            % evaluate Legendre polinomial for all points in the vector
            PVa=gsl_sf_legendre_Plm_e_mult(nn, m, cosThetaV);
            
            cosPhi = cos(m*PhiV); 
            sinPhi = sin(m*PhiV); 
            
            F = RV; 
            nF=1;
            nc=1;
            for n=nn
                for c=(nF+1):n, F = F .* RV; end
                nF=n;
                i=find(obj.nm(1,:)==n & obj.nm(2,:)==m);

%                 for j=1:3
%                     if abs(obj.AB(j,i))>eps && (j == 1 || j == 3)
%                         DV(j,:) = DV(j,:) + F.*PVa(nc,:).*cosPhi * obj.AB(j,i);
%                     end
%                     if abs(obj.AB(j,i))>eps && j == 2
%                         DV(j,:) = DV(j,:) + F.*PVa(nc,:).*sinPhi * obj.AB(j,i);
%                     end
%                 end
                for j=1:3
                    if abs(obj.A3(j,i))>eps
                        DV(j,:) = DV(j,:) + F.*PVa(nc,:).*cosPhi * obj.A3(j,i);
                    end
                    if abs(obj.B3(j,i))>eps
                        DV(j,:) = DV(j,:) + F.*PVa(nc,:).*sinPhi * obj.B3(j,i);
                    end
                end
                nc=nc+1;
            end
        end
        DV=obj.R0_mm*DV;
%         disp('DV');
%         disp(DV); %Verstehe ich nicht, wieso werrden alle moeglichen DV berechnet ... und am Ende Ueberschrieben???
    end
    
    
    function update_AB_from_sph_coeff(obj,x0)        
%         coeff =  obj.coeff; %old sph_coeff
%         obj.coeff.value  %change
%         coeff = x0
%         for i = 1:size(x0)
%             obj.coeff(i).value = x0(i);
%         end
        t = num2cell(x0);
        [obj.coeff.value] = t{:};
        coeff = obj.coeff;
%         nmax = max(nmax, mmax);
%         obj.coeff = coeff
        coeffDim = obj.coeffDim;
          
        %_initCoeff();
        Alpha_x = zeros(coeffDim);
        Alpha_y = zeros(coeffDim);
        Alpha_z = zeros(coeffDim);
        Beta_x  = zeros(coeffDim);
        Beta_y  = zeros(coeffDim);
        Beta_z  = zeros(coeffDim);

        %/**************************************************************************/
        %/****** organize coefficient values ******/
        for n = 1:size(x0,2)
            switch [coeff(n).A_or_B coeff(n).xyz]
                case 'Ax'
                    Alpha_x(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                case 'Ay'
                    Alpha_y(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                case 'Az'
                    Alpha_z(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                case 'Bx'
                    Beta_x(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                case 'By'
                    Beta_y(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                case 'Bz'
                    Beta_z(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                otherwise
                    error('ERROR: unrecognized coefficient string: ''%c%c''\n', coeff(n).A_or_B, coeff(n).xyz);
            end
        end
                
        
        % newer even more efficient structures optimized for vector calculations
        % nm : list of used pairs of n amd m indexes
        % A3 and B3 : all Alphas and Betas
        ii=find((abs(Alpha_x')>eps) | (abs(Beta_x')>eps) ...
               |(abs(Alpha_y')>eps) | (abs(Beta_y')>eps) ...
               |(abs(Alpha_z')>eps) | (abs(Beta_z')>eps));
        ni=length(ii);        
        obj.nm=zeros(2,ni);
        obj.A3=zeros(3,ni);
        obj.B3=zeros(3,ni);
        for i=1:ni
            [mp1,np1]=ind2sub([coeffDim coeffDim],ii(i));
            n=np1-1;
            m=mp1-1;
            obj.nm(1,i)=n;
            obj.nm(2,i)=m;
            
            % Siemens's normalization of spheical harmonics
            if m>0
                normfact = (-1)^m*sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m))); % TODO: checkme!
            else
                normfact=1;
            end
            
            obj.A3(1,i)=normfact*Alpha_x(np1,mp1);
            obj.B3(1,i)=normfact*Beta_x(np1,mp1);
            obj.A3(2,i)=normfact*Alpha_y(np1,mp1);
            obj.B3(2,i)=normfact*Beta_y(np1,mp1);
            obj.A3(3,i)=normfact*Alpha_z(np1,mp1);
            obj.B3(3,i)=normfact*Beta_z(np1,mp1);
            %test NW AB coeff Matrix
%             obj.AB(1,i)=normfact*Alpha_x(np1,mp1);
%             obj.AB(2,i)=normfact*Beta_y(np1,mp1);
%             obj.AB(3,i)=normfact*Alpha_z(np1,mp1);
            
        end
    end 
    
    function XDistorted = lab_to_mri(obj, XLinear)
        % transform coordinates from lab (mm) frame to MRI image frame
        % XL is the position vector in mm in the laboratory (scanner) frame
        % XD is the corresponding position in distorted MRI image coordinates in mm
        XDistorted = XLinear+spharm_evaluate(obj, XLinear); % NW?: fehlt da der erste term oder wieso haben wir die Addition am Anfang?
    end
    
    
    
    function XLinear = mri_to_lab(obj, XDistorted)
        % transform coordinates from image frame to lab frame
        % XD is the position vector in distorted MRI image coordinates in mm
        % XL is the corresponding position in mm in the laboratory (scanner) frame
        % XD may have two dimensions to store a list of vectors
        XLinear=zeros(size(XDistorted));
        for i=1:size(XDistorted,2)
            [XLinear(:,i),fval]=fminsearch(@(x) vecnorm(obj.lab_to_mri(x)-XDistorted(:,i)), XDistorted(:,i));
%             fval
        end
    end
    
    
%/*!
%\fn void GradUnwarp::create_transtable(VOL_GEOM *vg, MATRIX *vox2ras, MATRIX *inv_vox2ras)
%\brief This method creates GCAM (m3z transform table) for given VOL_GEOM using loaded gradient file.
%\param vg          - input VOL_GEOM struct
%\param vox2ras     - input vox2ras matrix
%\param inv_vox2ras - input inverse of vox2ras, ras2vox matrix
%*/
%void GradUnwarp::create_transtable(VOL_GEOM *vg, MATRIX *vox2ras, MATRIX *inv_vox2ras)
    function create_transtable(obj, vg)
        % create transformation table for the given voxel geometry (vg)
        % structure based on the spherical (solid) harmonic expansion
        % initialized/loaded in the class
        
        vox2ras = MRIxfmCRS2XYZ(vg); 
        
        %   // 1. create GCAM struct gcam, update gcam->image, gcam->atlas
        %   // 2. for each unwarped crs    % what is crs???
        %   //      calculate unwarped ras   %what is ras???
        %   //      convert unwarped ras to unwarped from RAS to LAI   what
        %   is LAI???
        %   //      call spharm_evaluate() to calculate delta ras in LAI
        %   //      warped ras = unwarped ras + delta ras
        %   //      convert warped ras from LAI to RAS
        %   //      calculate warped crs (fcs, frs, fss)
        %   //      update GCAM nodes - gcam->nodes are indexed by unwarped (c, r, s), 
        %   //                          (origx, origy, origz) is unwarped crs, 
        %   //                          (x, y, z) is warped crs
        
        m_RAS=ones(4,vg.len1,vg.len2,vg.len3);
        [m_RAS(1,:,:,:) m_RAS(2,:,:,:) m_RAS(3,:,:,:)]=ndgrid(0:(vg.len1-1),0:(vg.len2-1),0:(vg.len3-1));
        ras2lai=eye(4);
        ras2lai(1,1)=-1;
        ras2lai(3,3)=-1;
        m_RAS=ras2lai*vox2ras*m_RAS(:,:);
        
        % m_RAS1=permute(reshape(m_RAS,4,vg.len1,vg.len2,vg.len3),[2,3,4,1]);
        % figure; imagesc(m_RAS1(:,:,ceil(end/2),1)); title('voxel coordinate c1');colorbar;
        % figure; imagesc(m_RAS1(:,:,ceil(end/2),2)); title('voxel coordinate c2');colorbar;
        % figure; imagesc(m_RAS1(:,:,ceil(end/2),3)); title('voxel coordinate c3');colorbar;
        % squeeze(m_RAS1(end/2,end/2,ceil(end/2),:))
        
        DV = obj.spharm_evaluate(m_RAS(1:3,:));
        
        m_DistortedRAS = m_RAS;
        m_DistortedRAS(1:3,:)=m_DistortedRAS(1:3,:)+DV;
        m_DistortedCRS = inv(vox2ras) * ras2lai * m_DistortedRAS; 
        obj.gcam=permute(reshape(m_DistortedCRS(1:3,:),[3 vg.len1 vg.len2 vg.len3]),[2 3 4 1]); 
        
        % amplitude correction (determinant of the Jacobian of the transformation)
        nd=ndims(obj.gcam(:,:,:,1));
        u=cell(nd);
        if nd==3
            for i=1:3
                [u{i,1},u{i,2},u{i,3}] = gradient(obj.gcam(:,:,:,i));
            end
            obj.jd=zeros(size(u{1,1})); % unbelievable, I have to do the determinant calculation myself!
            for i=1:3
                for j=1:3
                    obj.jd = obj.jd + (-1)^(i+j)*u{i,j}.*(u{mod(i,3)+1,mod(j,3)+1}.*u{mod(i+1,3)+1,mod(j+1,3)+1} - ...
                                                          u{mod(i+1,3)+1,mod(j,3)+1}.*u{mod(i,3)+1,mod(j+1,3)+1});
                end
            end
            obj.jd=abs(obj.jd);
        else
            for i=1:2
                [u{i,1},u{i,2}] = gradient(obj.gcam(:,:,1,i));
            end
            obj.jd = abs(u{1,1}.*u{2,2} - ...
                         u{2,1}.*u{1,2});
        end
    end
    
% /*!
% \fn void GradUnwarp::read_siemens_coeff(const char *gradfilename)
% \brief This method reads and parses the gradient coefficient file.
% \param gradfilename   - input gradient coefficient file
% */
% void GradUnwarp::read_siemens_coeff(const char *gradfilename)
    function [coeff] = read_siemens_coeff(obj, gradfilename)
        %// check if gradfile has extension .grad

        %// open gradfilename
        fgrad = fopen(gradfilename, 'r');
        if fgrad < 0
            error('ERROR: could not read the coefficient file %s\n', gradfilename);
        end
        
        fprintf('==> reading coefficients from gradient coil file %s\n' , gradfilename) ;
        %_skipCoeffComment();
        
        [line, pos] = obj.skipComments(fgrad);
        line=strip(line);

        %// hard-coded limits:
        %// length of each entry    - 1024
        %char coeffline[1024];

        %// skip the next line. (It contains an information about the system type.)
        %fgets(coeffline, sizeof(coeffline), fgrad);

        fprintf('==> reading system type string from coeff.grad file...\n');
        fprintf('%s\n', line);

        % // check if first paramline contains "win_"; and, if so, parse
        % fgets(coeffline, sizeof(coeffline), fgrad);
        % if (strncmp(coeffline, " win_", 5) == 0)
        % {
        % 
        % // parse into the four parameters (these don't seem to be used anywhere...)
        % //[dum, iThreshLo, dum, iThreshUp, dum, iAlgoTyp, dum, iWinDummy] ...
        % //	= strread(paramline, " %10c%d%13c%d%13c%d%14c%d;");
        % 
        % // read next line
        % fgets(coeffline, sizeof(coeffline), fgrad);
        % }
        line=fgets(fgrad); % win_ ...
        line=fgets(fgrad);

        %// only extract radius and ignore rest
        R0=sscanf(line, '%f');
        obj.R0_mm = R0 * 1000;  %// transform R0 from m to mm
  
        % // read next line, which contains gradient system mode "(0 = typ. tunnel magnet system; 1 = typ. open magnet system)"
        % fgets(coeffline, sizeof(coeffline), fgrad);
        % int CoSyMode;
        line=fgets(fgrad);
        CoSyMode = sscanf(line, '%d');
        %fprintf("CoSyMode = %d\n", CoSyMode);

        %// skip the next 5 lines
        for ind = 1:5
            line=fgets(fgrad);
        end
        
        %// remember the starting of spherical harmonic coefficients
        fpos = ftell(fgrad);
        
        %// first pass: get the entry counts
        coeffCount=0;
        line=fgets(fgrad);
        while ischar(line)
            % int len = strlen(coeffline);
            % char* ptr = coeffline;
            % char* endptr = ptr + len;
            % 
            % // skip leading spaces and tabs, also empty lines
            % while (ptr != endptr && (*ptr == ' ' || *ptr == '\t' || *ptr == '\r' || *ptr == '\n'))
            %   ptr++;
            % 
            % if (*ptr == '\0')
            %   continue;
            line = strip(line);
            if isempty(line)
                line=fgets(fgrad);
                continue;
            end
                
            %if (getenv("GRADUNWARP_PRN_GRADCOEFF_READ"))
            %  printf("(first pass) entry #%d: %s\n", coeffCount, coeffline);
            
            coeffCount = coeffCount + 1;
            line=fgets(fgrad);
        end

        %coeff = new COEFF[coeffCount];
        coeff(1:coeffCount) = struct('num',0,'A_or_B','A','n',0,'m',0,'value',0.0,'xyz','x');
        nmax=0;
        mmax=0;

        %// rewind file pointer
        fseek(fgrad, fpos, -1);

        %/******************** second pass **************************/
        %/****** begin reading spherical harmonic coefficients ******/
        coeffCount=0;
        line=fgets(fgrad);
        while ischar(line)
            % int len = strlen(coeffline);
            % char* ptr = coeffline;
            % char* endptr = ptr + len;
            % 
            % // skip leading spaces and tabs, also empty lines
            % while (ptr != endptr && (*ptr == ' ' || *ptr == '\t' || *ptr == '\r' || *ptr == '\n'))
            %   ptr++;
            % 
            % if (*ptr == '\0')
            %   continue;
            line = strip(line);
            if isempty(line)
                line=fgets(fgrad);
                continue;
            end

            %if (getenv("GRADUNWARP_PRN_GRADCOEFF_READ"))
            %  printf("(second pass) entry #%d: %s\n", coeffCount, coeffline);

            A=sscanf(line, "%d %c(%d, %d) %f %c");
            if 6~=length(A)
                error('error parsing input file, line %s\n',line);
            end
            coeff(coeffCount+1).num=A(1);
            coeff(coeffCount+1).A_or_B=char(A(2));
            coeff(coeffCount+1).n=A(3);
            coeff(coeffCount+1).m=A(4);
            coeff(coeffCount+1).value=A(5);
            coeff(coeffCount+1).xyz=char(A(6));
            
            nmax = max(coeff(coeffCount+1).n, nmax);
            mmax = max(coeff(coeffCount+1).m, mmax);

            % if (getenv("GRADUNWARP_PRN_GRADCOEFF_READ"))
            %   printf("(second pass) %d %c (%d, %d) %f %c\n", coeff[coeffCount].num, 
            %          coeff[coeffCount].A_or_B, coeff[coeffCount].n, coeff[coeffCount].m, coeff[coeffCount].value, coeff[coeffCount].xyz);
            
            coeffCount = coeffCount + 1;
            line=fgets(fgrad);
        end

        fclose(fgrad);    

        nmax = max(nmax, mmax);
        obj.coeff = coeff
        coeffDim = nmax+1;
        obj.coeffDim = nmax+1;

        %_initCoeff();
        Alpha_x = zeros(coeffDim);
        Alpha_y = zeros(coeffDim);
        Alpha_z = zeros(coeffDim);
        Beta_x  = zeros(coeffDim);
        Beta_y  = zeros(coeffDim);
        Beta_z  = zeros(coeffDim);

        %/**************************************************************************/
        %/****** organize coefficient values ******/
        for n = 1:coeffCount
            switch [coeff(n).A_or_B coeff(n).xyz]
                case 'Ax'
                    Alpha_x(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                case 'Ay'
                    Alpha_y(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                case 'Az'
                    Alpha_z(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                case 'Bx'
                    Beta_x(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                case 'By'
                    Beta_y(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                case 'Bz'
                    Beta_z(coeff(n).n+1,coeff(n).m+1)=coeff(n).value;
                otherwise
                    error('ERROR: unrecognized coefficient string: ''%c%c''\n', coeff(n).A_or_B, coeff(n).xyz);
            end
        end
                
        % % new more efficient structures that avoid unnecessary calculations
        % ii_x=find((abs(Alpha_x(:))>eps) | (abs(Beta_x(:))>eps));
        % ni_x=length(ii_x);        
        % obj.abnm_x = struct('Alpha',0,'Beta',0,'n',0,'m',0);
        % for i=1:ni_x
        %     [n,m]=ind2sub([obj.coeffDim obj.coeffDim],ii_x(i));            
        %     obj.abnm_x(i).n=n-1;
        %     obj.abnm_x(i).m=m-1;
        %     obj.abnm_x(i).Alpha=obj.Alpha_x(n,m);
        %     obj.abnm_x(i).Beta=obj.Beta_x(n,m);
        % end
        % ii_y=find((abs(obj.Alpha_y(:))>eps) | (abs(obj.Beta_y(:))>eps));
        % ni_y=length(ii_y);        
        % obj.abnm_y = struct('Alpha',0,'Beta',0,'n',0,'m',0);
        % for i=1:ni_y
        %     [n,m]=ind2sub([obj.coeffDim obj.coeffDim],ii_y(i));            
        %     obj.abnm_y(i).n=n-1;
        %     obj.abnm_y(i).m=m-1;
        %     obj.abnm_y(i).Alpha=obj.Alpha_y(n,m);
        %     obj.abnm_y(i).Beta=obj.Beta_y(n,m);
        % end
        % ii_z=find((abs(obj.Alpha_z(:))>eps) | (abs(obj.Beta_z(:))>eps));
        % ni_z=length(ii_z);        
        % obj.abnm_z = struct('Alpha',0,'Beta',0,'n',0,'m',0);
        % for i=1:ni_z
        %     [n,m]=ind2sub([obj.coeffDim obj.coeffDim],ii_z(i));            
        %     obj.abnm_z(i).n=n-1;
        %     obj.abnm_z(i).m=m-1;
        %     obj.abnm_z(i).Alpha=obj.Alpha_z(n,m);
        %     obj.abnm_z(i).Beta=obj.Beta_z(n,m);
        % end
        
        % newer even more efficient structures optimized for vector calculations
        % nm : list of used pairs of n amd m indexes
        % A3 and B3 : all Alphas and Betas
        ii=find((abs(Alpha_x')>eps) | (abs(Beta_x')>eps) ...
               |(abs(Alpha_y')>eps) | (abs(Beta_y')>eps) ...
               |(abs(Alpha_z')>eps) | (abs(Beta_z')>eps));
        ni=length(ii);        
        obj.nm=zeros(2,ni);
        obj.A3=zeros(3,ni);
        obj.B3=zeros(3,ni);
        for i=1:ni
            [mp1,np1]=ind2sub([coeffDim coeffDim],ii(i));
            n=np1-1;
            m=mp1-1;
            obj.nm(1,i)=n;
            obj.nm(2,i)=m;
            
            % Siemens's normalization of spheical harmonics
            if m>0
                normfact = (-1)^m*sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m))); % TODO: checkme!
            else
                normfact=1;
            end
            
            obj.A3(1,i)=normfact*Alpha_x(np1,mp1);
            obj.B3(1,i)=normfact*Beta_x(np1,mp1);
            obj.A3(2,i)=normfact*Alpha_y(np1,mp1);
            obj.B3(2,i)=normfact*Beta_y(np1,mp1);
            obj.A3(3,i)=normfact*Alpha_z(np1,mp1);
            obj.B3(3,i)=normfact*Beta_z(np1,mp1);
            %test NW AB coeff Matrix
            obj.AB(1,i)=normfact*Alpha_x(np1,mp1);
            obj.AB(2,i)=normfact*Beta_y(np1,mp1);
            obj.AB(3,i)=normfact*Alpha_z(np1,mp1);
            
        end

        obj.AB_initialized = true;
    end

    function [nextLine,filePosBefore] = skipComments(obj,fid)
        %skipComments Read lines of skipping blank lines and comments.
        %   line=skipComments(fid) Read lines from valid file identifer and
        %   return the next non-comment line.
        
        filePosBefore=ftell(fid);
        line = fgetl(fid);
        while ischar(line) && (isempty(line) || line(1) == '#')
            filePosBefore=ftell(fid);
            line = fgetl(fid);
        end
        if ischar(line)
            nextLine = line;
        else
            filePosBefore = -1;
            nextLine = -1;
        end
    end
    
end
end

% local functon used in the class
function m=MRIxfmCRS2XYZ(mri)
% find transform from image coordinates 
% (column row slice, CRS) to XYZ
  m=eye(4);

  %/* direction cosine between rows scaled by
  %   distance between rows */
  m(1:3,1)=mri.v1_ras * mri.vsize1;

  %/* direction cosine between columns scaled by
  %   distance between colums */
  m(1:3,2)=mri.v2_ras * mri.vsize2;

  %/* direction cosine between slices scaled by
  %   distance between slices */
  m(1:3,3)=mri.v3_ras * mri.vsize3;

  %/* At this point, m = Mdc * D */

  %/* Col, Row, Slice at the Center of the Volume */
  Pcrs = ones(4,1);
  Pcrs(1,1)=mri.len1/2;
  Pcrs(2,1)=mri.len2/2;
  if mri.len3>1
    Pcrs(3,1)=mri.len3/2;
  else
    Pcrs(3,1)=0;
  end

  %/* XYZ offset the first Col, Row, and Slice from Center */
  PxyzOffset = m * Pcrs;

  %/* XYZ at the Center of the Volume is mri->c_r, c_a, c_s  */

  %/* The location of the center of the voxel at CRS = (0,0,0)*/
  m(1:3,4)=mri.c_ras' - PxyzOffset(1:3);
  
end

