##### connect to a node
```
srun --pty bash
```

##### see current jobs
```
squeue -u username
```

##### nb threads
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
du -sh
```

##### cancel all jobs
```
squeue -u username | awk '{print $1}' | tail -n+2 | xargs scancel
```

##### copy files to cluster
```
scp -r username@genologin.toulouse.inra.fr:/work/username/your_folder/ .
```

##### run Rmd scripts  
```
module load system/R-4.2.3_Miniconda3
Rscript -e "rmarkdown::render('example.Rmd')
```