clear all
close all

%% main function
% path = '/home/localwehkamp/distortion_matrix/distortion_matrix_test/';
% path = '/home/localwehkamp/distortion_matrix/distortion_matrix_test3_230721/';
path = '/home/localwehkamp/distortion_matrix/144/';

pose1_points = get_points(append(path, 'calibration_1.calib'));
%Write function that reads all calibration positions and returns 3 dim
%array
% number_of_calib_files = 38;
% number_of_calib_files = 24;
number_of_calib_files = 120;

point_clouds = pose1_points;
for i = 2:number_of_calib_files
    number = int2str(i);
    pose_points = get_points(append(path, 'calibration_', number ,'.calib'));
    point_clouds = cat(3,point_clouds,pose_points);
end

disp('points');
disp(pose1_points(:,1));
disp('points.shape');
disp(size(pose1_points));
disp('points.dtype');
disp(class(pose1_points));
%disp(rotation0.as_euler('zxy', degrees=True));
%disp(rotation0.as_matrix());

reference_points = point_clouds(:,:,1); %pose5_points  !!!change to first pose in future Version
points = point_clouds(:,:,:);
% points = point_clouds(:,:,10:30);
example_measurement = 12;

%%Initialize gradient unwarping class
gradUnwarp = GradUnwarpV(); 

%Read spherical harmonic coefficients for (gradient) distortion matrix 
[coeff] = gradUnwarp.read_siemens_coeff('coeff_AS82.grad');

x = [coeff.value];
[sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff(points(:,:,example_measurement), reference_points, gradUnwarp, x);

%% Creating figure
fig = figure('Units', 'pixels', 'Position', [0 0 1000 700]);
ax = axes('Parent', fig, 'DataAspectRatio', [1 1 1], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'Projection', 'perspective');
hold(ax, 'on');

% Creating plot
scatter3(ax, reference_points(:,1), reference_points(:,2), reference_points(:,3), 'g', 'x');
scatter3(ax, translated_rotated_points(:,1), translated_rotated_points(:,2), translated_rotated_points(:,3), 'r', 'x');
% title('Field

%% Optimize sph_coeffs that minimize sum_vecnorm_difference_pose

x0 = [coeff.value];
%%
% gradUnwarp.update_AB_from_sph_coeff(x0);
% [sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff(points, reference_points, gradUnwarp, x0);
measured_points = points(:,:,:);

f = @(x) find_sph_coeff(measured_points, reference_points, gradUnwarp, x);    %numeric function with vector input
options = optimset('TolFun', 1e-64,'PlotFcns',@optimplotfval);
[x,fval,exitflag,output] = fminsearch(f,x0,options);
x
%%
[sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff(points, reference_points, gradUnwarp, x);
%%
% %difference_pose = reference_points[valid_points] - translated_rotated_points[valid_points];
% mean_difference_pose = mean(difference_pose, 1);
% disp('mean_difference_pose');
% disp(mean_difference_pose);
% max_difference_pose = max(vecnorm(difference_pose'));
% disp('max_difference_pose_vecnorm');
% disp(max_difference_pose);


% [valid_points_1, invalid_points] = extract_valid_points_based_on_reference_distances(points, reference_points);
% disp('valid_points');
% disp(valid_points_1);
% disp('valid_points.shape');
% disp(size(valid_points_1));
% [translated_rotated_points, rotation, translation1, translation2] = match_points_to_origin(points, reference_points, valid_points_1); %(points, origin_points, iterations)

% % [valid_points, invalid_points] = extract_valid_points_after_transform_to_lab_coordinate(points, reference_points);
% %translated_points = points - translation1;
% %rotated_points = rotation.apply(translated_points);
% %translated_rotated_points = rotated_points + translation2;
% difference_pose_1 = reference_points - translated_rotated_points;
% disp('difference_pose_1');
% disp(difference_pose_1);
% %difference_pose = reference_points[valid_points] - translated_rotated_points[valid_points];
% mean_difference_pose = mean(difference_pose_1, 1);
% disp('mean_difference_pose');
% disp(mean_difference_pose);
% max_difference_pose = max(vecnorm(difference_pose_1'));
% disp('max_difference_pose_vecnorm');
% disp(max_difference_pose);



%% Difference in linear (Lab) space


% [valid_points, invalid_points] = extract_valid_points_based_on_reference_distances(points, reference_points);
% disp('valid_points');
% disp(valid_points);
% disp('valid_points.shape');
% disp(size(valid_points));
% [translated_rotated_points, rotation, translation1, translation2] = match_points_to_origin(points, reference_points, valid_points); %(points, origin_points, iterations)
% 
% % [valid_points, invalid_points] = extract_valid_points_after_transform_to_lab_coordinate(points, reference_points);
% %translated_points = points - translation1;
% %rotated_points = rotation.apply(translated_points);
% %translated_rotated_points = rotated_points + translation2;
% difference_pose = reference_points - translated_rotated_points;
% disp('difference_pose');
% disp(difference_pose);
% %difference_pose = reference_points[valid_points] - translated_rotated_points[valid_points];
% mean_difference_pose = mean(difference_pose, 1);
% disp('mean_difference_pose');
% disp(mean_difference_pose);
% max_difference_pose = max(vecnorm(difference_pose'));
% disp('max_difference_pose_vecnorm');
% disp(max_difference_pose);




%% bi ier

% field_offset_from_point_displacement(points, reference_points);
% field_offset_from_point_displacement(translated_rotated_points, reference_points);

% Creating figure
fig = figure('Units', 'pixels', 'Position', [0 0 1000 700]);
ax = axes('Parent', fig, 'DataAspectRatio', [1 1 1], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'Projection', 'perspective');
hold(ax, 'on');

% Creating plot
scatter3(ax, reference_points(:,1), reference_points(:,2), reference_points(:,3), 'g', 'x');
scatter3(ax, translated_rotated_points(:,1,example_measurement), translated_rotated_points(:,2,example_measurement), translated_rotated_points(:,3,example_measurement), 'r', 'x');
% title('Field

%%
% define function to extract 16 points form .scan file
% the points are in the order of the list, start after the keyword "positions"
function mat_list_points = get_points(file_name)
    file = fopen(file_name, 'r');
    line = fgetl(file);
    str_start = strfind(line,'positions')+13;
    str_end = strfind(line,'basisID')-5;
    str_points = line(str_start:str_end);
    list_points = strsplit(str_points, '],[');
    for i = 1:numel(list_points)
        split_list_points{i} = str2double(strsplit(list_points{i},','));
    end
    mat_list_points = cell2mat(split_list_points);
    mat_list_points = reshape(mat_list_points,[3,16]);
    mat_list_points = transpose(mat_list_points);
%     np_points = cell2mat(split_list_points);
end


function [rotation] = kabsch_algorithm(A,B) %align_vectors %Thank you Marco
    centroid_A = mean(A, 2);
    centroid_B = mean(B, 2);
    A_centered = A; %- centroid_A;% !!! Achtung habe ich jetzt 2 mal centriert??? ANscheinend ja
    B_centered = B; %- centroid_B;
    covariance = A_centered * B_centered';
    [U, ~, V] = svd(covariance);
    rotation = V * U';
% B_aligned = rotation * B;
end

% function [translated_rotated_point_clouds] = match_point_clouds_to_origin(point_clouds, iso_center_points, valid_points_in_clouds) %valid_points_in_clouds??
%     for i in 3:%?shape size of point_clouds
%         points = point_clouds(i)
%         translated_rotated_point_clouds[i] = match_points_to_origin(points, iso_center_points, valid_points);
% end

% function [translated_rotated_points, rotation, mean_center_new, mean_center_origin] = match_points_to_origin(points, iso_center_points, valid_points)
function [translated_rotated_points, rotation, mean_center_new, mean_center_origin] = match_points_to_origin(points, iso_center_points, valid_points)
    mean_center_origin = mean(iso_center_points(valid_points,:));
    origin_points = iso_center_points(valid_points,:) - mean_center_origin;
    mean_center_new = mean(points(valid_points,:));
    translated_points = points(valid_points,:) - mean_center_new;
    rotation = kabsch_algorithm(origin_points', translated_points');
    translated_points = points - mean_center_new;
    rotated_points = rotation * translated_points';
    translated_rotated_points = rotated_points' + mean_center_origin;
end

function [translated_rotated_points, rotation, mean_center_new, mean_center_origin] = match_points_to_origin_old(points, iso_center_points, valid_points)
    mean_center_origin = mean(iso_center_points(valid_points,:));
    origin_points = iso_center_points(valid_points,:) - mean_center_origin;
    mean_center_new = mean(points(valid_points,:));
    translated_points = points(valid_points,:) - mean_center_new;
    rotation = kabsch_algorithm(origin_points', translated_points');
    translated_points = points - mean_center_new;
    rotated_points = rotation * translated_points';
    translated_rotated_points = rotated_points' + mean_center_origin;
end

function [valid_points, invalid_points] = extract_valid_points_based_on_threshold_distances(points, reference_points,threshold) %Validation based on small position distortion assumption in center region 
    pdist_ref = pdist(reference_points, 'euclidean');
    pdist_new = pdist(points, 'euclidean'); %speed up with pdist(X,"fasteuclidean",CacheSize=10);
%     disp('pdist_new.size', pdist_new(:,:,1))
    pdist_diff = pdist_ref - pdist_new;
    sf = squareform(pdist_diff); %????
    distance_threshold = threshold; %0.003;
    [sf_index,~] = find(sf < distance_threshold); % threshold to exclude points based on changes in the relative distances
    occurences = histcounts(sf_index);
    valid_points = find(occurences > 3);
%     size(valid_points)
    if size(valid_points,2) <= 2
        error('Warning, less or equal than 2 points excced distance threshold!');
    end
    invalid_points = find(occurences <= 3);
end

% function [valid_points, invalid_points] = extract_valid_points_after_transform_to_lab_coordinate(points, reference_points)
%     points=points'*1000;
%     % xxr=reshape(gradUnwarp.mri_to_lab(xx(:,1)),size(xx(:,1)));
%     xxr=reshape(gradUnwarp.mri_to_lab(points),size(points))
%     
%     valid_points = 0;
%     invalid_points = 0;
% end
    
% function [test] = get_spatial_deviation_from_linear_gradient(points,invalid_points,transformation,rotation...)
%     
%     test = points;
% end

% function field_offset_from_point_displacement(translated_rotated_points, reference_points)
%     difference_pose = reference_points - translated_rotated_points;
%     dist = sqrt(sum(difference_pose.^2, 2));
%     B0_freq = 123e6; %42.58 % in [MHz]
%     field_offset = dist * B0_freq; % in [Hz]
%     x_field_offset = difference_pose(:,1) * B0_freq; % in [Hz]
%     y_field_offset = difference_pose(:,2) * B0_freq; % in [Hz]
%     z_field_offset = difference_pose(:,3) * B0_freq; % in [Hz]
% end

% function sph_harm()
%     sph_harm(m, n, theta, phi, 'vectorized');
% end
% function [sum_vecnorm_difference_pose, translated_rotated_points] = old_find_sph_coeff(points, reference_points, gradUnwarp, x0)
%     gradUnwarp.update_AB_from_sph_coeff(x0);
% 
%     %Transform to linear "lab" space using initial sph_coeffs (wrong, but best initial guess)
%     xxd = points'*1000; % why multiply 1000 ???
%     xxl=reshape(gradUnwarp.mri_to_lab(xxd),size(xxd));
%     lab_points = xxl'/1000;
%     % lab_points = points
% 
% %     difference_pose = reference_points - lab_points;
% %     disp('difference_pose');
% %     disp(difference_pose);
% 
%     % valid_points_1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
%     threshold = 0.00001;
%     threshold = 0.001;  % should I increase the threashold with each itteration??? Unnecessay
%     [valid_points, invalid_points] = extract_valid_points_based_on_threshold_distances(lab_points, reference_points,threshold);
% %     disp('valid_points');
% %     disp(valid_points);
% %     disp('valid_points.shape');
% %     disp(size(valid_points));
%     [translated_rotated_points, rotation, translation1, translation2] = match_points_to_origin(lab_points, reference_points, valid_points); %(points, origin_points, iterations)
% 
%     % translated_rotated_points = translated_rotated_points + 0.001
%     %translated_points = points - translation1;
%     %rotated_points = rotation.apply(translated_points);
%     %translated_rotated_points = rotated_points + translation2;
%     difference_pose = reference_points - translated_rotated_points;
% %     disp('difference_pose');
% %     disp(difference_pose);
% 
%     sum_vecnorm_difference_pose = sum(vecnorm(difference_pose'));
% 
% %     disp('sum_vecnorm_difference_pose');  
% %     disp(sum_vecnorm_difference_pose);
%     
% %     fig = figure('Units', 'pixels', 'Position', [0 0 1000 700]);
% %     ax = axes('Parent', fig, 'DataAspectRatio', [1 1 1], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'Projection', 'perspective');
% %     hold(ax, 'on');
% % 
% %     % Creating plot
% %     scatter3(ax, pose5_points(:,1), pose5_points(:,2), pose5_points(:,3), 'g', 'x');
% %     scatter3(ax, translated_rotated_points(:,1), translated_rotated_points(:,2), translated_rotated_points(:,3), 'r', 'x');
% end

% function [sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff(points, reference_points, gradUnwarp, x0)
%     %% Update AB from spherical harmonic coefficients listed in "x0"
%     gradUnwarp.update_AB_from_sph_coeff(x0);
% 
%     %% Transform to linear "lab" space using initial sph_coeffs (wrong, but best initial guess)
%     xxd = points'*1000; % why multiply 1000 ???
%     xxl=reshape(gradUnwarp.mri_to_lab(xxd),size(xxd));
%     lab_points = xxl'/1000;
%     % lab_points = points
% 
%     % valid_points_1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
%     threshold = 0.001;  % should I increase the threashold with each itteration??? Unnecessay?
%     [valid_points, invalid_points] = extract_valid_points_based_on_threshold_distances(lab_points, reference_points,threshold); %give an error when less then 3 points are valid
% %     disp('valid_points');
% %     disp(valid_points);
% %     disp('valid_points.shape');
% %     disp(size(valid_points));
% %     [translated_rotated_points, rotation, translation1, translation2] = match_points_to_origin(lab_points, reference_points, valid_points); %(points, origin_points, iterations)
%     [translated_rotated_points] = match_points_to_origin(lab_points, reference_points, valid_points); %(points, origin_points, iterations)
% %     [translated_rotated_point_clouds] = match_point_clouds_to_origin(point_clouds, iso_center_points, valid_points_in_clouds)
% %     for i = 1:numel(list_points)%oder size shape...
% %         difference_pose_clouds(i) = reference_points - translated_rotated_point_clouds(i);
% %     end
% 
%     difference_pose = reference_points - translated_rotated_points;
% %     disp('difference_pose');
% %     disp(difference_pose);
% 
% %     for i = 1:numel(list_points)%oder size shape...
% %         sum_vecnorm_difference_pose = sum(vecnorm(difference_pose'));
% %     end
%     sum_vecnorm_difference_pose = sum(vecnorm(difference_pose'));
% 
% 
% end


% function [sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff_old2(points, reference_points, gradUnwarp, x0)
%     %% Update AB from spherical harmonic coefficients listed in "x0"
%     gradUnwarp.update_AB_from_sph_coeff(x0);
% 
%     %% Transform to linear "lab" space using initial sph_coeffs (wrong, but best initial guess)
%     number_of_point_clouds = numel(points(1,1,:));
% %     xxd = points(:,:,28)'*1000; % why multiply 1000 ???
% %     xxd = pagetranspose(points)*1000; % why multiply 1000 ???
% %     xxl = reshape(gradUnwarp.mri_to_lab(xxd),size(xxd));
%     
%     xxl = zeros(size(points));
%     lab_points = zeros(size(points));
%     for i = 1:number_of_point_clouds
%         xxd = points(:,:,i)'*1000; % why multiply 1000 ???
% %         xxl(:,:,i)=reshape(gradUnwarp.mri_to_lab(xxd),size(xxd));
%         xxl = reshape(gradUnwarp.mri_to_lab(xxd),size(xxd));
%         lab_points(:,:,i) = xxl'/1000;
%     end
% %     lab_points = xxl'/1000;
% %     lab_points = pagetranspose(xxl)/1000;
%     % lab_points = points
% 
%     % valid_points_1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
%     threshold = 0.001;  % should I increase the threashold with each itteration??? Unnecessay?
% %     [valid_points, invalid_points] = extract_valid_points_based_on_threshold_distances(lab_points, reference_points,threshold); %give an error when less then 3 points are valid
% 
% %     valid_points = zeros()
% %     invalid_points = zeros()
%     A = [];
%     for i = 1:number_of_point_clouds %Attention this does not work with reasonable threshold as it leads to a change in matrix size, possible solution boolean matrix..
%         [valid_points, invalid_points] = extract_valid_points_based_on_threshold_distances(lab_points(:,:,i), reference_points,threshold);
%         A = cat(3,A,valid_points); 
%     end
% %     size(A);
%         %     disp('valid_points');
% %     disp(valid_points);
% %     disp('valid_points.shape');
% %     disp(size(valid_points));
% %     [translated_rotated_points, rotation, translation1, translation2] = match_points_to_origin(lab_points, reference_points, valid_points); %(points, origin_points, iterations)
% %     [translated_rotated_point_clouds] = match_point_clouds_to_origin(point_clouds, iso_center_points, valid_points_in_clouds)
% 
% %     [translated_rotated_points] = match_points_to_origin(lab_points, reference_points, valid_points); %(points, origin_points, iterations)
%     translated_rotated_points = zeros(size(points));
%     for i = 1:number_of_point_clouds
%         [translated_rotated_points(:,:,i)] = match_points_to_origin(lab_points(:,:,i), reference_points, valid_points); %(points, origin_points, iterations)
%     end
%     
%     
% %     for i = 1:numel(list_points)%oder size shape...
% %         difference_pose_clouds(i) = reference_points - translated_rotated_point_clouds(i);
% %     end
%     
%     difference_pose = reference_points - translated_rotated_points;
% %     disp('difference_pose');
% %     disp(difference_pose);
% 
%     sum_vecnorm_difference_pose = 0;
%     for i = 1:number_of_point_clouds %oder size shape...
%         sum_vecnorm_difference_pose = sum_vecnorm_difference_pose + sum(vecnorm(difference_pose(:,:,i)'));
%     end
% %     sum_vecnorm_difference_pose = sum(vecnorm(difference_pose'));
% 
% 
% end
function [sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff(measured_points, reference_points, gradUnwarp, x0)
    %% Update AB from spherical harmonic coefficients listed in "x0"
    gradUnwarp.update_AB_from_sph_coeff(x0);

    %% Transform to linear "lab" space using initial sph_coeffs (wrong, but best initial guess)
    number_of_point_clouds = numel(measured_points(1,1,:));
    
    xxl = zeros(size(measured_points));
    lab_points = zeros(size(measured_points));
%     tic
    parfor i = 1:number_of_point_clouds
        xxd = measured_points(:,:,i)'*1000; % why multiply 1000 ???
        xxl = reshape(gradUnwarp.mri_to_lab(xxd),size(xxd));
        lab_points(:,:,i) = xxl'/1000;
    end
%     toc
    threshold = 0.1;  % should I increase the threashold with each itteration??? Unnecessay?

%     valid_points = zeros()
%     invalid_points = zeros()
    A = [];
    for i = 1:number_of_point_clouds %Attention this does not work with reasonable threshold as it leads to a change in matrix size, possible solution boolean matrix..
        [valid_points, invalid_points] = extract_valid_points_based_on_threshold_distances(lab_points(:,:,i), reference_points,threshold);
        A = cat(3,A,valid_points); 
    end

    translated_rotated_points = zeros(size(measured_points));
    parfor i = 1:number_of_point_clouds
        [translated_rotated_points(:,:,i)] = match_points_to_origin(lab_points(:,:,i), reference_points, valid_points); %(points, origin_points, iterations)
    end
    
    difference_pose = reference_points - translated_rotated_points;

%     sum_vecnorm_difference_pose = 0;
%     sums_vecnorm_difference_pose = zeros(number_of_point_clouds);
    for i = 1:number_of_point_clouds
        sums_vecnorm_difference_pose(i) = sum(vecnorm(difference_pose(:,:,i)'));
    end
    sum_vecnorm_difference_pose = sum(sums_vecnorm_difference_pose)
    
%     sum_vecnorm_difference_pose = 0;
%     for i = 1:number_of_point_clouds
%         sum_vecnorm_difference_pose = sum_vecnorm_difference_pose + (sum(vecnorm(difference_pose(:,:,i)')))^2;
%     end
%     parfor i = 1:number_of_point_clouds
%         sum_vecnorm_difference_pose(i) = sum_vecnorm_difference_pose + (sum(vecnorm(difference_pose(:,:,i)')))^2;
%     end

end