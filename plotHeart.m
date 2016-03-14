obj = wf_load('Human Heart lowerhalf.obj');
[ vox, res ] = voxelize(obj);
plot3D(vox, res);