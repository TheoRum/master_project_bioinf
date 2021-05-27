#!/bin/bash

for f in $HOME/hd/SLE/DATA_7T/rsfmri_200909/Nifti_master/sub-*; do
sub=$(echo $f | awk -F 'sub-' '{print $2}' | cut -d"/" -f2)

echo subject $sub is started ...

docker run -itd --name ppfMRIprep -v $HOME/hd/SLE/DATA_7T/rsfmri_200909/Nifti_master:/data:ro -v $HOME/hd/SLE/DATA_7T/rsfmri_200909/ppfMRIprep/derivatives:/out -v $HOME/hd/SLE/DATA_7T/rsfmri_200909/ppfMRIprep/workdir:/worktemp -v $FREESURFER_HOME:/FS_folder poldracklab/fmriprep:20.2.0 /data /out/out participant --participant_label $sub --nproc 24 --bold2t1w-dof 9 --force-bbr --output-spaces T1w MNI152NLin6Asym --fs-no-reconall --use-aroma --stop-on-first-crash --fs-license-file /FS_folder/license.txt -w /worktemp --write-graph -v

docker wait ppfMRIprep
echo ...subject $sub is finished.
docker logs ppfMRIprep > $HOME/hd/SLE/DATA_7T/rsfmri_200909/$sub.log
docker container rm ppfMRIprep
echo docker removed
done
