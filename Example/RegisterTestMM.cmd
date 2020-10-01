::Adjust ANN path
set ANN_PATH=E:\Computations\Reconstruction\HighRes\Data\Program\ann_bin
set PATH=%PATH%;%ANN_PATH%

::Adjust program path
call ..\VS2010\MM_Restricted\Release\MM_Restricted.exe ..\Model\All_30_7.rmm ..\Model\MeanFace.off  ..\Model\All_Lmks.txt ..\Example\stereo_pointcloud.off ..\Example\stereo_pointcloud_landmarks.txt ..\Example\stereo_pointcloud_fitting.off

pause

