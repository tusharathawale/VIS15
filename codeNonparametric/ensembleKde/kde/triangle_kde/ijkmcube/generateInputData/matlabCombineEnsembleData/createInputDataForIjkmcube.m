% inputdir:  Directory containing all ensemble members, e.g., fuel_input
% source:  Names of files in input directory have form i_source.vud. Here
% source is 'noisy'
% numMembers: number of ensemble members
% outputdir:  Directory to store mean, i.e, target.nrrd, and all members
% with kde bandwidth 2_target.nrrd.
% target: Name of NRRD file, e.g., if target string is 'fuel', the name of
% the nrrd file will be fuel.nrrd.
% sample function call: createInputDataForIjkmcube('fuel_input/', 'noisy', 10, 'fuel_output/', 'fuel')
function createInputDataForIjkmcube(inputdir, source,  numMembers, outputdir, target)
addpath('nrrd_read_write_rensonnet/');
v = loadVUD(strcat(inputdir,'0_', source,'.vud'));
[d1,d2,d3] = size(v);

% Collect all members into a single array
ensemble = zeros(d1, d2, d3, numMembers+1);
ensemble(:,:,:,1) = v;
for i = 1:numMembers-1
     ensemble(:,:,:,i+1) = loadVUD(strcat(inputdir,num2str(i),'_',source,'.vud'));
end

meanVol = zeros(d1, d2, d3);

% Compute bandwidth of Kernel density estimation and store it
for i=1:d1
    for j=1:d2
        for k=1:d3
             x = squeeze(ensemble(i,j,k,1:numMembers));
            
             
             % The bandwidth correction for input kernel choice, e.g.,
             % uniform, triangle, Epanechnikov is applied in ijkmcube (C++) code
             % by referring to Table 1 of Appendix
             % standard deviation
             std_dev = std(x);
             % interquartile range
             r = iqr(x);
             % data spread
             spread = min (std_dev, (1/1.34)*r); 
             meanVol(i,j,k) = mean(x);  
             %Bandwidth using the Silverman's rule of thumb with Gaussian kernel (Refer to Table 1 and Section 2 of Appendix)
			 bandwidth = ((1.0)^-0.4) * ((0.5/sqrt(pi))^0.2) * ((3.0^8.0 * (22.0/7.0)^-0.5)^-0.2) * (numMembers^-0.2) * spread;
             ensemble(i,j,k,numMembers+1) = bandwidth;
             
        end
    end
end


ensemble = permute(ensemble,[4 1 2 3]);

% Store all ensemble members
% Struct specifying header of NRRD file
s1 = struct;
s1.data = ensemble;
s1.sizes = [11, 64, 64, 64];
s1.encoding = 'raw';
s1.dimension = 4;
s1.type = 'double';

headerInfo = nhdr_nrrd_write(strcat(outputdir,'2_',target,'.nrrd'), s1, true);

% Store mean
% Struct specifying header of NRRD file
s2 = struct;
s2.data = meanVol;
s2.sizes = [64, 64, 64];
s2.encoding = 'raw';
s2.dimension = 3;
s2.type = 'double';

headerInfo = nhdr_nrrd_write(strcat(outputdir, target,'.nrrd'), s2, true);


