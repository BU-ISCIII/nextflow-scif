# Nextflow tutorial workflow using SCIF and Singularity/Docker
This repository implements the Nextflow tutorial workflow using the Scientific Filesystem (SCIF) and Docker/Singularity to provide a reproducible research environment. Is an adaptation of [snakemake.scif](https://github.com/sci-f/snakemake.scif), but using nextflow instead. Also I have just discovered [this try](https://github.com/vsoch/rnatoy.scif/edit/master/README.md) of @vsoch and I am adapting her cool documentation to this case.

These are the steps to follow in order to develop a new reproducible and transparent workflow:

- **SCIF recipes**: **Discoverability** and **Transparency** by installing our software in the container via a [Scientific Filesystem](https://sci-f.github.io). SCIF also lets us install the same dependencies across container technologies.
- **Singularity/Docker container**: this will handle the **Reproducibility** of our pipeline and the software it is going to need. Moreover Docker and Singularity HUB makes deployment and sharing a lot easier (**TO DO**).
- **Nextflow**enables scalable and reproducible scientific workflows managing software containers natively like a charm. It allows the adaptation of pipelines written in the most common scripting languages. Moreover its fluent DSL simplifies the implementation and the deployment of complex parallel and reactive workflows on clouds and clusters, being able to abstract the users both from *where* and even *how* they compute without any concern about the multiple dependencies a scientific workflow can have.
- **Testing**: circleCI (**TO DO**)

So...Let's begin..

## 1. Recipes
When building pipelines, you can think of it like baking a cake. We have entire recipes for creating our final products (containers), and within those recipes ingredients (software) that we need to add. In this first part, we will talk about the three recipes in this repository, the [Dockerfile](Dockerfile) for the Docker container, the [Singularity](Singularity) recipe for a Singularity container, and the [app_recipes/*.scif](app_recipes) recipes. 

### The Scientific Filesystem Recipe
A scientific filesystem is useful because it allows me to write one recipe for my various software, and then install easily in different containers or on my host. How do you know when you find a recipe? When you find a recipe for a scientific filesystem (SCIF), you will see a file with extension *.scif. For example, in this repository:

 - [app_recipes/*.scif](app_recipes) are the recipes for the scientific filesystem that will be installed in both the Docker and Singularity containers to be run with nextflow. SCIF is flexible in that there can be **many** different internal applications defined in one file, however if we want we can put them in individual files and install them equivalently, as we show in this example. 
 
 For example, given the apps "samtools", "bwa" and "bcftools" and using a three recipes files , I would install like:

```
scif install app_recipes/samtools_v1.9_centos7.scif
scif install app_recipes/bwa_v0.7.17_centos7.scif
scif install app_recipes/bcftools_v1.9_centos7.scif
```

but I could also define the different applications in only one file and install them the next way.

```
scif install rnatoy.scif
```

This level of modularity is up to the user. Programatically, this install would be equivalent. This is a purely "human information" decision. If I am building a single container to share with my paper, I would opt to only need one file. If I am providing a builder service to users and want to easily install recipes in a modular fashion, I would want to do the second, which is this case. Here we are preparing a model for developing new workflows as quickly as possible, so app_recipes should contain a collection of app recipes that can be combined in the Singularity/Docker recipe for creating the needed container for the workflow being developed. In fact, probably we will convert app_recipes to a Github repository, so it will be used as remote repository in all the workflows we will develop, making its update easier.

This means that, for any SCIF recipe and a host of interest (Docker or Singularity container, or your computer) you can install the same recipes. 

## The Docker Recipe
[Dockerfile](Dockerfile) is the recipe for building our container. You will also notice the installation is simple - we start with a container base that was equivalently used by the creator of the pipeline with system / host dependencies, and then simply install the SCIF recipe to it. That comes down to these three commands:

```
RUN pip install scif                           # Install scif from pypi
ADD app_recipes/*.scif /opt                    # Add the recipe to the container
RUN scif install /opt/rnatoy.scif               # Install it to the container
```

## The Singularity Recipe
The recipe file for a Singularity container is the file [Singularity](Singularity). The format of the recipe file is different, but installing the scientific filesystem, again from the recipes [app_recipes](app_recipes) is performed with the same commands:

```
pip install scif                                           # Install scif from pypi
scif install /opt/samtools_v1.9_centos7.scif               # Install it to the container
```

The only missing command to add the recipe to the container is because Singularity recipes allow you to do this in a `%files` section that functions like a cp <source> <dest> (if not dest supplied the file will be copied in ```/```)

```
%files
    samtools_v1.9_centos7.scif /opt
```
