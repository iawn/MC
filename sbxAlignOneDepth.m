function sbxAlignOneDepth(outputBase,useTheseFiles,hi,wi,MD,depth)
%wrapper function for a single depth sbx align
%created by ian, jan 2017

OPB=[outputBase '_depth_' num2str(depth)];
disp('Beginning MC')
sbxalignmastermulti_3D(useTheseFiles,1,[OPB '.align'],hi{depth},wi{depth},MD,depth);

save([OPB '.align'],'outputBase','-append');
save([OPB '.align'],'file','-append');
save([OPB '.align'],'loadStart','-append');
save([OPB '.align'],'numToLoad','-append');
save([OPB '.align'],'hi','wi','-append');
try
save([OPB '.align'],'rect','-append');
end
%sbxcomputeci(useTheseFiles,[outputBase '.align'],hi{depth},wi{depth}); %Takes about 10 minutes, eats up a ton of RAM

end