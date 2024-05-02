clear all
close all

%% main function
path = '../144/';

pose1_points = get_points(append(path, 'calibration_1.calib'));
%Write function that reads all calibration positions and returns 3 dim
%array
number_of_calib_files = 122;
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
%%
reference_points = point_clouds(:,:,1); %pose5_points  !!!change to first pose in future Version
%%
p1column = point_clouds(:,:,1:11);
p2column = point_clouds(:,:,14:24);
p3column = point_clouds(:,:,25:34);
p4column = point_clouds(:,:,39:48);
p5column = point_clouds(:,:,49:58);
p6column = point_clouds(:,:,64:72);
p7column = point_clouds(:,:,73:81);
p8column = point_clouds(:,:,89:96);
p9column = point_clouds(:,:,97:103);
p10column = point_clouds(:,:,116:120);
p11column = point_clouds(:,:,121:122);
measured_points = cat(3,p1column,p2column,p3column,p4column,p5column,p6column,p7column,p8column,p9column,p10column,p11column);

example_measurement = 15;

%%Initialize gradient unwarping class
gradUnwarp = GradUnwarpV(); 

%Read initial spherical harmonic coefficients for (gradient) distortion matrix 
x = struct2array(load("x0","-mat"));

[sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff(measured_points(:,:,example_measurement), reference_points, gradUnwarp, x);
% [sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff_multi(measured_points, reference_points, gradUnwarp, x);

% %% bi ier
% xxdr = reference_points'*1000;
% xxlr = reshape(gradUnwarp.mri_to_lab(xxdr), size(xxdr));
% reference_points = xxlr'/1000;

%% Creating figure
fig = figure('Units', 'pixels', 'Position', [0 0 1000 700]);
ax = axes('Parent', fig, 'DataAspectRatio', [1 1 1], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'Projection', 'perspective');
hold(ax, 'on');

% Creating plot
scatter3(ax, reference_points(:,1), reference_points(:,2), reference_points(:,3), 'g', 'x');
scatter3(ax, translated_rotated_points(:,1), translated_rotated_points(:,2), translated_rotated_points(:,3), 'r', 'x');
% scatter3(ax, translated_rotated_points(:,1,example_measurement), translated_rotated_points(:,2,example_measurement), translated_rotated_points(:,3,example_measurement), 'r', 'x');

% title('Field

%% Optimize sph_coeffs that minimize sum_vecnorm_difference_pose

x0 = struct2array(load("x0","-mat"));
% x0 = struct2array(load("x0_new11","-mat"));
% x0 = x0*0
%%
% gradUnwarp.update_AB_from_sph_coeff(x0);
% [sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff(points, reference_points, gradUnwarp, x0);


% f = @(x) find_sph_coeff(measured_points(:,:,example_measurement), reference_points, gradUnwarp, x);    %numeric function with vector input
f = @(x) find_sph_coeff_multi(measured_points, reference_points, gradUnwarp, x);    %numeric function with vector input
% options = optimset('TolFun', 1e-64,'PlotFcns',@optimplotfval);
options = optimset('MaxFunEvals',1e12,'MaxIter',1e12,'TolFun', 1e-64,'TolX', 1e-64,'PlotFcns',@optimplotfval);
[x,fval,exitflag,output] = fminsearch(f,x0,options);
x
save('x0_thresh0_00005.mat','x');
%%
[sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff(measured_points(:,:,example_measurement), reference_points, gradUnwarp, x);
% [sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff_multi(measured_points, reference_points, gradUnwarp, x);



%% Creating figure
fig = figure('Units', 'pixels', 'Position', [0 0 1000 700]);
ax = axes('Parent', fig, 'DataAspectRatio', [1 1 1], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'Projection', 'perspective');
hold(ax, 'on');

% Creating plot
scatter3(ax, reference_points(:,1), reference_points(:,2), reference_points(:,3), 'g', 'x');
scatter3(ax, translated_rotated_points(:,1), translated_rotated_points(:,2), translated_rotated_points(:,3), 'r', 'x');
% scatter3(ax, translated_rotated_points(:,1,example_measurement), translated_rotated_points(:,2,example_measurement), translated_rotated_points(:,3,example_measurement), 'r', 'x');

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


function [rotation] = kabsch_algorithm(A,B) %align_vectors
%     centroid_A = mean(A, 2);
%     centroid_B = mean(B, 2);
    A_centered = A; %- centroid_A;% !!! Achtung habe ich jetzt 2 mal centriert??? ANscheinend ja
    B_centered = B; %- centroid_B;
    covariance = A_centered * B_centered';
    [U, ~, V] = svd(covariance);
    rotation = V * U';
% B_aligned = rotation * B;
end


function [translated_rotated_points, rotation, mean_center_new, mean_center_origin] = match_points_to_origin(points, iso_center_points, valid_points)
    valid_points = logical(valid_points);
    mean_center_origin = mean(iso_center_points(valid_points,:));
    origin_points = iso_center_points(valid_points,:) - mean_center_origin;
    mean_center_new = mean(points(valid_points,:));
    translated_points = points(valid_points,:) - mean_center_new;
    rotation = kabsch_algorithm(origin_points', translated_points');
    translated_points = points - mean_center_new;
    rotated_points = rotation * translated_points';
    translated_rotated_points = rotated_points' + mean_center_origin;
end

function [valid_points] = extract_valid_points_based_on_proximity_to_iso(points)
    distance_to_iso = vecnorm(points,2,2);
    [B,I] = mink(distance_to_iso,3);
    valid_points = zeros(size(distance_to_iso));
    valid_points(I) = true;
end

function [valid_points, invalid_points] = extract_valid_points_based_on_threshold_distances(points, reference_points,distance_threshold) %Validation based on small position distortion assumption in center region 
    pdist_ref = pdist(reference_points, 'euclidean');
    pdist_new = pdist(points, 'euclidean'); %speed up with pdist(X,"fasteuclidean",CacheSize=10);
    pdist_diff = pdist_ref - pdist_new;
    sf = squareform(pdist_diff); % changes in distances of each point to the other points in the two point sets 
    [sf_index,~] = find(abs(sf) <= distance_threshold); % threshold to exclude points based on changes in the relative distances
    occurences = histcounts(sf_index);
    valid_points = occurences >= 3;
    if sum(valid_points) <= 2
        disp('Warning, less or equal than 2 points excced distance threshold!');
        distance_threshold = 0.0005;
        [sf_index,~] = find(abs(sf) <= distance_threshold); % threshold to exclude points based on changes in the relative distances
        occurences = histcounts(sf_index);
        valid_points = occurences >= 3;
        if sum(valid_points) <= 2
            valid_points = extract_valid_points_based_on_proximity_to_iso(points);
        end
    end
end

function [valid_points, invalid_points] = extract_3_valid_points(points, reference_points) %Validation based on small position distortion assumption in center region 
    pdist_ref = pdist(reference_points, 'euclidean');
    pdist_new = pdist(points, 'euclidean'); %speed up with pdist(X,"fasteuclidean",CacheSize=10);
    pdist_diff = pdist_ref - pdist_new;
    sf = squareform(pdist_diff); % changes in distances of each point to the other points in the two point sets 
    [sf_index,~] = find(abs(sf) <= distance_threshold); % threshold to exclude points based on changes in the relative distances
    occurences = histcounts(sf_index);
    valid_points = occurences >= 3;
    if sum(valid_points) <= 2
        disp('Warning, less or equal than 2 points excced distance threshold!');
    end
end


function [sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff_multi(measured_points, reference_points_d, gradUnwarp, x0)
    %% Update AB from spherical harmonic coefficients listed in "x0"
    gradUnwarp.update_AB_from_sph_coeff(x0);
    x0 %just to see whats going on

    %% Transform to linear "lab" space using initial sph_coeffs (wrong, but best initial guess)
    number_of_point_clouds = numel(measured_points(1,1,:));
    
    xxdr = reference_points_d'*1000;
    xxlr = reshape(gradUnwarp.mri_to_lab(xxdr), size(xxdr));
    reference_points = xxlr'/1000;

    xxl = zeros(size(measured_points));
    lab_points = zeros(size(measured_points));
    tic
    parfor i = 1:number_of_point_clouds
        xxd = measured_points(:,:,i)'*1000; % why multiply 1000 ???
        xxl = reshape(gradUnwarp.mri_to_lab(xxd),size(xxd));
        lab_points(:,:,i) = xxl'/1000;
    end
    toc
%     threshold = 0.000097;  %0.002start ; 0.0002 should I increase the threashold with each itteration???
    % threshold = 0.00005; 
    threshold = 0.0005;

    valid_points = zeros(1,16,number_of_point_clouds);
    for i = 1:number_of_point_clouds %Attention this does not work with reasonable threshold as it leads to a change in matrix size, possible solution boolean matrix..
        [valid_points(:,:,i)] = extract_valid_points_based_on_threshold_distances(lab_points(:,:,i), reference_points,threshold);
%         [valid_points(:,:,i)] = extract_valid_points_based_on_proximity_to_iso(lab_points(:,:,i));
    end

    translated_rotated_points = zeros(size(measured_points));
%     parfor i = 1:number_of_point_clouds
    for i = 1:number_of_point_clouds
        [translated_rotated_points(:,:,i)] = match_points_to_origin(lab_points(:,:,i), reference_points, valid_points(:,:,i)); %(points, origin_points, iterations)
    end
    
    difference_pose = reference_points - translated_rotated_points;

%     sum_vecnorm_difference_pose = 0;
%     sums_vecnorm_difference_pose = zeros(number_of_point_clouds);
    for i = 1:number_of_point_clouds
        sums_vecnorm_difference_pose(i) = sum(vecnorm(difference_pose(:,:,i)'));
    end
    sum_vecnorm_difference_pose = sum(sums_vecnorm_difference_pose)
    

end

function [sum_vecnorm_difference_pose, translated_rotated_points] = find_sph_coeff(points, reference_points_d, gradUnwarp, x0)
    %% Update AB from spherical harmonic coefficients listed in "x0"
    gradUnwarp.update_AB_from_sph_coeff(x0);

    %% Transform to linear "lab" space using initial sph_coeffs (wrong, but best initial guess)
    xxdr = reference_points_d'*1000;
    xxlr = reshape(gradUnwarp.mri_to_lab(xxdr), size(xxdr));
    reference_points = xxlr'/1000;

    xxd = points'*1000; % why multiply 1000 ???
    xxl=reshape(gradUnwarp.mri_to_lab(xxd),size(xxd));
    lab_points = xxl'/1000;
    % toc
    % lab_points = points

    threshold = 0.00015;  % should I decrease the threashold with each itteration??? 
    [valid_points] = extract_valid_points_based_on_threshold_distances(lab_points, reference_points,threshold); %give an error when less then 3 points are valid
%     [valid_points] = extract_valid_points_based_on_proximity_to_iso(lab_points);
    %     [valid_points, invalid_points] = extract_valid_points_based_on_threshold_distances(lab_points, reference_points,threshold); %give an error when less then 3 points are valid

    %     disp('valid_points');
%     disp(valid_points);
%     disp('sum(valid_points)');
%     disp(sum(valid_points));
%     [translated_rotated_points, rotation, translation1, translation2] = match_points_to_origin(lab_points, reference_points, valid_points); %(points, origin_points, iterations)
    [translated_rotated_points] = match_points_to_origin(lab_points, reference_points, valid_points); %(points, origin_points, iterations)
%     [translated_rotated_point_clouds] = match_point_clouds_to_origin(point_clouds, iso_center_points, valid_points_in_clouds)
%     for i = 1:numel(list_points)%oder size shape...
%         difference_pose_clouds(i) = reference_points - translated_rotated_point_clouds(i);
%     end

    difference_pose = reference_points - translated_rotated_points;
%     disp('difference_pose');
%     disp(difference_pose);

    sum_vecnorm_difference_pose = sum(vecnorm(difference_pose'));


end