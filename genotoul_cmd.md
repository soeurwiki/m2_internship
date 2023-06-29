##### connect to a node
```
srun --pty bash
```

##### see current jobs
```
squeue -u satia
```

##### nb coeurs
```
cat /proc/cpuinfo | grep processor | wc -l 
```

##### list empty files
```
find -type f -empty
```

##### nb of files
```
ls | wc -l
```

##### size of directory
```
du -sh work/
```

##### cancel all jobs
```
squeue -u satia | awk '{print $1}' | tail -n+2 | xargs scancel
```

##### copy files to cluster
```
scp -r satia@genologin.toulouse.inra.fr:/work/satia/folder/ .
```