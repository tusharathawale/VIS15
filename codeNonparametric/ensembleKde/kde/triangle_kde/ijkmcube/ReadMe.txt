% ----------------------------------------------------------------------------------
%   Title: Isosurface Visualization of Data with Nonparametric Models for Uncertainty
%   Authors: Tushar M. Athawale, Elham Sakhaee, and Alireza Entezari	
%   Date: 06/26/2019
% ---------------------------------------------------------------------------------------

Overview:

The code is provided for marching cubes algorithm in ensemble representing uncertain data for isosurface extraction using nonparametric noise models. The provided code assumes the triangle kernel for nonparametric density estimation from uncertain data. The code may be modified for uniform/Epanechnikov kernel assumption (see directions below). The given code works on ensemble data, and it assumes the vertex-based classification for topology prediction (Sec 3.1). The geometry extraction is similar to one explained in Sec. 3.2. 

The code for non-local means version (Sec 4), edge-based topological classification (Sec 3.1), and experiments for single scalar field data (Fig. 11, 12) are not provided, but will be made available soon..

———————————————————————————————————————————
For ijkmcube and geomview installation instructions:

Refer to the document ijkmcube_installation.txt

Steps for generating visualizations:

1) Preprocess uncertain data to feed it as input to the ijkmcube source code (C++):

        - Navigate to codeNonparametric/ensembleKde/kde/triangle_kde/ijkmcube/generateInputData/matlabCombineEnsembleData/ in terminal

        - Place all ensemble members (VUD format)in the input directory, such that their names have format i_filename.vud with i = 0 .. numOfMembers-1. For example, we placed all ensemble members for the fuel dataset into a directory fuel_input with their names as i_noisy.vud.

        - Run the provided MATLAB script createInputDataForIjkmcube.m using command createInputDataForIjkmcube('fuel_input/', 'noisy', 10, 'fuel_output/', 'fuel')

        - Running of MATLAB script will generate two files: fuel.nrrd and 2_fuel.nrrd in the output directory, i.e., fuel_output. The fuel.nrrd contains mean of ensemble members, and 2_fuel.nrrd contains all members stored together with bandwidth for kernel density estimation assuming a Gaussian kernel (the bandwidth correction is applied later in the ijkmcube code if a kernel other than the Gaussian, such as, uniform, triangle, or Epanechnikov, is used for kernel density estimation).
          
Thus, if a*b*c is a grid resolution input data with k ensemble members, fuel.nrrd has a resolution a*b*c and 2_fuel.nrrd has a resolutions a*b*c*(k+1) [+1 is for storing kernel density estimation bandwidth]

        - copy two files, i.e., fuel.nrrd and 2_fuel.nrrd, in the directory codeNonparametric/ensembleKde/kde/triangle_kde/ijkmcube/
        
2) We have provided the modified version of the IJKMCUBE (http://web.cse.ohio-state.edu/research/graphics/isotable/documentation/ijkmcube.html) code is provided for isosurface extraction in uncertain data assuming nonparametric distribution models.

        - Again, for instructions on IJKMCUBE installation, please refer to the document ijkmcube_installation.txt. The installation should generate an executable with name “ijkmcube”.
       
        - Now run one of the following commands. Note that the kernel bandwidth correction for triangle kernel is applied in the file ijkispoly.txx inside a function intersects_poly (Refer to Appendix to learn mathematical formulae for applying the kernel bandwidth correction. Also see the code tree provided at bottom of this page.).   
        a) ./ijkmcube -interpolate "alpha_uncertainty" 32 fuel.nrrd  (resolve ambiguous cases of the marching cubes in uncertain data using the mean of uncertain data)
        b) ./ijkmcube -interpolate "alpha_uncertainty" -topology "adecider" 32 fuel.nrrd (resolve ambiguous cases of the marching cubes in uncertain data using the probabilistic midpoint decider (Sec. 3.3))

3) The provided code considers the triangle base kernel for nonparametric model (Although, we always used the uniform base kernel for resolving the marching cubes face ambiguity cases.). If you would like to use the kernel other than the triangle, e.g., uniform or Epanechnikov for nonparametric density estimation, you will need to make the following three changes.

        - open ijkispoly.txx, and update bandwidth correction formulae for the needed kernel (comment out line for the triangle kernel, and uncomment line for the kernel of your choice) inside function intersects_poly (See below detailed code tree to get idea of functionality of intersects_poly function)

        - open ijkinterpolate.txx, use handle corresponding to the kernel of your choice (uni, tri, or ep) inside the function kde_alpha_uncertainty_coord

        - open ijkinterpolate.txx, update bandwidth correction for the needed kernel (comment out line for the triangle kernel, and uncomment line for the kernel of your choice)
      

—————————————————————————————————————————————————————————————
Ijkmcube code tree for this project:

IJKMCUBE (http://web.cse.ohio-state.edu/research/graphics/isotable/documentation/ijkmcube.html):

We first modified the IJKMCUBE code to handle uncertain data and topological uncertainty (vertex-based classification, Sec. 3.1) in the Marching cubes algorithm.
We integrated IJKMCUBE with our code file functions.cxx to resolve topologically ambiguous cases in the marching cubes algorithm. Finally, we integrated the IJKMCUBE with our code files, namely, kdeWithUniformKernel, kdeWithTriangleKernel, and kdeWithEpanechnikov kernel to derive the geometric uncertainty in recontructed isosurfaces for uncertain data (Sec. 3.2).

ijkmcube_main.cxx (main function)

  ijkmcube.cxx (call to the marching_cubes function for certain/uncertain input data)  
	 
     1) Determine the most probable isosurface topology in uncertain data (Sec. 3.1) and resolve ambiguous topological cases.  

       a) The path (a) is taken when the probabilistic midpoint decider (Sec. 3.3) is not used to resolve ambiguous topology cases of the marching cubes in uncertain data or when no ' -topology "adecider" ' option is used with the ijkmcube command (see step 2 above)
       
	-> ijkmcube_extract.cxx (Determine isosurface topology within each grid cube by computing grid vertex signs with vertex-based classification approach) (Function: extract_iso_simplices, extract_iso_simplices_binary) 

             -> ijkispoly.txx (code for determining which cubes of scalar field grid are crossed by an isosurface followed by the application of vertex-based classification for isosurface-crossed cubes (Section 3.1 in the paper)) (Functions names: extract_isopatch_from_mesh_poly, intersects_poly (to check if cell is crossed), compute_isotable_index (vertex-based classification))

             -> This route will resolve ambiguous cases of the marching cubes algorithm assuming the mean of uncertain data

      OR

      b) The path (b) is taken when probabilistic midpoint decider (Sec. 3.3) is used to resolve ambiguous topology cases of the marching cubes in uncertain data or when ' -topology "adecider" ' option is used with the ijkmcube command (see step 2 above)
     
	 ->  ijkmcube_extract.cxx (Determine isosurface topology within each grid cube by computing grid vertex signs with vertex-based classification approach. If the most probable configuration represents the ambiguous configuration of the marching cubes algorithm, the ambiguity is resolved using the probabilistic midpoint decider (Sec. 3.3 in the paper). Note that the cell body ambiguities in the current code are still resolved using the mean of uncertain data, but the face ambiguities are resolved using the probabilistic midpoint decider.) (Functions: extract_iso_simplices, extract_iso_simplices_adecider, extract_from_cube_adecider, extract_from_ambig_cube_adecider, compute_signs_of_facet_deciders)
                  
	-> functions.cxx (Perform convolution of nonparametric distributions with the uniform kernel as the base) (Function: MidPointWrapper)

    2) Extract the isosurface geometry for decided topology and the positional variance in the geometry (Sec. 3.2)
       
      -> ijkmcube_sub.cxx (Functions: merge_and_position_vertices, position_iso_vertices_linearB_alpha, position_iso_vertices_linear_binaryB_alpha)
              
          -> ijkinterpolate.txx (Function: kde_alpha_uncertainty_coord)       
         

      

 