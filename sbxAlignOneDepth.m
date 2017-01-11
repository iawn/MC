function sbxAlignOneDepth(outputBase,useTheseFiles,hi,wi,MD,depth)
%wrapper function for a single depth sbx align
%created by ian, jan 2017

OPB=[outputBase '_depth_' num2str(depth)];
disp('Beginning MC')
sbxalignmastermulti_3D(useTheseFiles,1,[OPB '.align'],hi{depth},wi{depth},MD,depth);
%return %%debug code

end