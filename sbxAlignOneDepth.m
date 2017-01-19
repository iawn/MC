function sbxAlignOneDepth(outputBase,useTheseFiles,hi,wi,MD,depth,useRed)
%wrapper function for a single depth sbx align
%created by ian, jan 2017

computeci = 1; %should be turned on for nomal use, turned off to make function go faster for debugging

OPB=[outputBase '_depth_' num2str(depth)];
disp('Beginning MC')
sbxalignmastermulti_3D(useTheseFiles,computeci,[OPB '.align'],hi{depth},wi{depth},MD,depth,useRed);
%return %%debug code

end