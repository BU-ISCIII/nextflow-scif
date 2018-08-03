# Nextflow tutorial workflow using SCIF and Singularity/Docker
This repository implements the Nextflow tutorial workflow using the Scientific Filesystem (SCIF) and Docker/Singularity to provide a reproducible research environment. Is an adaptation of [snakemake.scif](https://github.com/sci-f/snakemake.scif), but using nextflow instead.

**NOTE:** I don't have configured the installation for Docker, I feel a little lazy right now.  This is a proof of concept and I have decided to use the approximation in master.

This is an alternative configuration with some pros and cons:
- In this configuration we have installed nextflow inside the container. This means that have remove one dependency, we can get this pipeline to work only with Singularity or Docker installed.
- Also we have separated apps environments, apps binaries are not added to the path in the container recipe, but in nextflow app, if you use any app independently you will only add that app env variables. However when you run nextflow app all the app bin folders will be added to the PATH.
- However I found a main con for this approximation, nextflow needs to have write permissions in the folder it is being executed, which in this case is /scif/app. I had to cd to /scif/data prior to nextflow execution and the directory with your data (your $PWD) must be bind to /scif/data in your singularity command.This makes the execution command look like this:

```
singularity run --bind $PWD:/scif/data ./nexflow-scif-alt.simg run nextflow run main.nf
```
which I personally don't like so much. 

Also, you could do it pointing to Github repository:

```
singularity run --bind $PWD:/scif/data ./nexflow-scif-alt.simg run nextflow run -r feature/nextflow-scif-alt BU-ISCIII/nextflow-scif
```
It is not a bad solution either, but I think it is more difficult to understand, and I can't explain why but I dislike not using singularity and docker as nextflow natively intends to xD.

If you want to see the configuration I prefer please check master branch.
