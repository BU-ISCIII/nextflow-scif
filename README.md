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
scif install nextflow-scif.scif
```

This level of modularity is up to the user. Programatically, this install would be equivalent. This is a purely "human information" decision. If I am building a single container to share with my paper, I would opt to only need one file. If I am providing a builder service to users and want to easily install recipes in a modular fashion, I would want to do the second, which is this case. Here we are preparing a model for developing new workflows as quickly as possible, so app_recipes should contain a collection of app recipes that can be combined in the Singularity/Docker recipe for creating the needed container for the workflow being developed. In fact, probably we will convert app_recipes to a Github repository, so it will be used as remote repository in all the workflows we will develop, making its update easier.

This means that, for any SCIF recipe and a host of interest (Docker or Singularity container, or your computer) you can install the same recipes. 

SCIF makes the installation to be so tidy, it makes me wanna cry, a scif folder is created in ```/``` of your container and the tree directory will look like this:
```
/scif/apps/
     samtools/
        bin/
        lib/
        scif/
            runscript.help
            runscript
            env/
                01-base.sh
                90-environment.sh
     bwa/
     ....
```

Also SCIF will provide a separate variable environment for each app (even we are going to mess it a little bit in order to get nextflow functioning), and some cool ENV variables $SCIF_APP(BIN|ENV|...) that seem promising for a lot of functionality. You can read more on [Singularity docs](https://www.sylabs.io/guides/2.5.1/user-guide/reproducible_scif_apps.html) and [SCIF docs](https://sci-f.github.io/).


### The Docker Recipe
[Dockerfile](Dockerfile) is the recipe for building our container. You will also notice the installation is simple - we start with a container base that was equivalently used by the creator of the pipeline with system / host dependencies. In this case, I have used a base docker container of centos 7, this is because we use this SO in our lab and our HPC, and it is useful for me to have homogeneous installation recipes in both containers, my workstation and HPC (at least for now that we are figuring all this stuff out).

Then simply install the SCIF recipe to it. That comes down to these three commands:

```
RUN pip install scif                           # Install scif from pypi
ADD app_recipes/*.scif /opt                    # Add the recipe to the container
RUN scif install /opt/rnatoy.scif               # Install it to the container
```
Now a little about the "mess", you can see that we added export PATH variable with bin folder for each app, this should not be necessary because app/bin folder is added to the PATH automaticaly when the app is run, **BUT** nextflow functionality expects to have access to the executables in the path or you have to provide the full path.
```
# Docker Recipe
ENV PATH=${PATH}:/scif/apps/samtools/bin 
# Singularity Recipe
echo 'export PATH=${PATH}:/scif/apps/samtools/bin' >> $SINGULARITY_ENVIRONMENT
```
In order to take advantage of complete nextflow functionality I decided to assume this mistake in the use of SCIF, which will make the apps envs not being totally separated. You can read the [Notes](#notes) section for the pros and cons of this approximation and an insight of another posibility I don't like so much (or I consider less useful for our situation).

**NOTE**: it would be cool to use $SCIF_APPBIN_samtools, instead of the full path ```/scif/apps/samtools/bin```, but I tried and this variable must be set after PATH env variable, because I can't get it working.

### The Singularity Recipe
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

## 2. Containers build and usage
Once we have our recipe or in order to create our installation commands with ease, we need to know how to build and interact with containers.

### Shell to a container
You can interact with a container in many different ways and for differente purposes. For example, during the development of the container, I took a strategy to start with a base, interactively shell into it, and test installation and running of things. To do this we can directly use a centos container from Dockerhub, singularity is what I prefer:

```bash
sudo singularity shell -w docker://centos:latest
```

Or docker if you want:
```bash
docker run -it centos
```

### Build a container
Once you have your recipe prepared you need to build your image:

With Singularity
```bash
sudo singularity build nextflow-scif.img Singularity
```

Or with Docker:
```bash
## . points to the path where docker will import the context for image creation. If you are in the repository path you have only to do:
docker build -t nextflow-scif .
# Dockerfile and context will be automatically imported to the docker daemon.
```
### Running a container
Command run, in Singularity and Docker will execute environment variables and the script configured in ```%runscript``` section in [Singularity](Singularity) or ENTRYPOINT in [Dockerfile](Dockerfile). In this case only scif is configured, being the command that will let us interact with the applications installed in the container.
**NOTE:** We can't define an ENTRYPOINT in docker because Nextflow doesn't support docker containers with entrypoint at the moment, it seems in its [docs](https://www.nextflow.io/docs/latest/docker.html#docker-page) that have implemented some twist but is already deprecated. This means that for interacting with the applications we have to explicitly use scif command after docker run.

```bash
# Docker
docker run nextflow-scif scif

# Singularity
./nextflow-scif
singularity run nextflow-scif
```

The output will be something like that

```
Scientific Filesystem [v0.0.71]
usage: scif [-h] [--debug] [--quiet] [--writable]
            
            {version,pyshell,shell,preview,help,install,inspect,run,apps,dump,exec}
            ...

scientific filesystem tools

optional arguments:
  -h, --help            show this help message and exit
  --debug               use verbose logging to debug.
  --quiet               suppress print output
  --writable, -w        for relevant commands, if writable SCIF is needed

actions:
  actions for Scientific Filesystem

  {version,pyshell,shell,preview,help,install,inspect,run,apps,dump,exec}
                        scif actions
    version             show software version
    pyshell             Interactive python shell to scientific filesystem
    shell               shell to interact with scientific filesystem
    preview             preview changes to a filesytem
    help                look at help for an app, if it exists.
    install             install a recipe on the filesystem
    inspect             inspect an attribute for a scif installation
    run                 entrypoint to run a scientific filesystem
    apps                list apps installed
    dump                dump recipe
    exec                execute a command to a scientific filesystem
```

#### Inspecting Applications
The strength of SCIF is that it will always show you the applications installed in a container, and then provide predictable commands for inspecting, running, or otherwise interacting with them. For example, if I find the container, without any prior knowledge I can reveal the applications inside:

```
# Docker
docker run nextflow-scif scif apps
# Singularity
./nextflow-scif apps
# or
singularity run nextflow-scif apps
```
The output will look like this:
```
  samtools
  bwa
  bcftools
```

We can also obtain help about a particular application

```bash
# Docker
docker run nextflow-scif scif help samtools
# Singularity
./nextflow-scif help samtools
# or
singularity run nextflow-scif help samtools
```
The output will look like this:
```
    Samtools
```

and then inspecting

```
# Docker
docker run nextflow-scif scif inspect samtools
# Singularity
./nextflow-scif inspect samtools
# or
singularity run nextflow-scif inspect samtools
```
Output:
```
{
    "samtools": {
        "apprun": [
            "    exec samtools \"$@\""
        ],
        "apphelp": [
            "   Samtools"
        ],
        "applabels": [
            "VERSION 1.9",
        ]
    }
}
```

The creator of the container didn't write any complicated scripts to have this happen - the help text is just a chunk of text in a block of the recipe. The labels that are parsed to json, are also just written easily on two lines. This means that the creator can spend less time worry about exposing this. If you can write a text file, you can make your applications programatically parseable. You can check how this is done adding free text to ```%apphelp``` and ```%applabel``` in [samtools_v1.9_centos7.scif](app_recipes/samtools_v1.9_centos7.scif).

**NOTE:** ```%appenv``` can be also used for environment variables particular of each application. We are not using this because of nextflow functionality. See [Notes](#notes).

#### Interacting with Applications
I can easily shell into the container in the context of an application, meaning that the
environment is sourced, etc. **THIS IS NOT WORKING: USING SHELL OR EXEC ON A APP RUNS %RUNSCRIPT OF THE APP, TRY REMOVING THE EXPORT PATHS**

```
# Docker
docker run -it nextflow-scif scif shell samtools
# Singularity
./nextflow-scif shell samtools
```

Output:

```
[samtools] executing /bin/bash 
root@d002e338b88b:/scif/apps/samtools# env | grep PATH
LD_LIBRARY_PATH=/scif/apps/samtools/lib
PATH=/scif/apps/samtools/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
```

Notice how I'm in the app's context (in it's application folder) and that it's bin is added to the path? I can also shell in without a specific application context, but still have all the SCIF [global variables](https://sci-f.github.io/spec-v1#environment-namespace) available to me.

### Using the container without the apps

We can still shell into the "general" container:

```
# Docker
docker run -it nextflow-scif shell
# Singularity
./nextflow-scif shell
# or
singularity shell nextflow-scif

```
Output:
```
WARNING No app selected, will run default ['/bin/bash']
executing /bin/bash 
root@055a34619d17:/scif# ls
apps
data
```
Or thanks our "mess" in which we have exported bin app folders to the path we can execute each app directly from singularity container:
```
### This are the commands that nextflow will run internally.
# Docker 
docker run nextflow-scif samtools
# Singularity
singularity exec nextflow-scif samtools
```

#### Running Applications
Before we get into creating a pipeline, look how easy it is to run an application. Without scif, we would have to have known that samtools is installed, and then executed the command to the container. But with the scientific filesystem, we discovered the app (shown above) and then we can just run it. The `run` command maps to the entrypoint, as was defined by the creator:

```
# Docker
docker run nextflow-scif run samtools
# Singularity 
./nextflow-scif run samtools
```
Output:

```
Program: samtools (Tools for alignments in the SAM format)
Version: 0.1.18 (r982:295)

Usage:   samtools <command> [options]

Command: view        SAM<->BAM conversion
         sort        sort alignment file
         mpileup     multi-way pileup
         depth       compute the depth
         faidx       index/extract FASTA
         tview       text alignment viewer
         index       index alignment
         idxstats    BAM index stats (r595 or later)
         fixmate     fix mate information
         flagstat    simple stats
         calmd       recalculate MD/NM tags and '=' bases
         merge       merge sorted alignments
         rmdup       remove PCR duplicates
         reheader    replace BAM header
         cat         concatenate BAMs
         targetcut   cut fosmid regions (for fosmid pool only)
         phase       phase heterozygotes

[samtools] executing /bin/bash /scif/apps/samtools/scif/runscript
```

And executing any command in the context of the application is possible too: **THIS IS NOT WORKING**.
This is equal to the issue shown with shell, due to out little twist for get nextflow working, here we have all the apps added to the PATH, but this should be the correct functioning of SCIF:

```
# Docker
docker run nextflow-scif scif exec samtools env | grep PATH
# Singularity
./nextflow-scif exec samtools env | grep PATH
```
Output:
```
LD_LIBRARY_PATH=/scif/apps/samtools/lib
PATH=/scif/apps/samtools/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
```

Note that for above, you will get more output with the Singularity container, as it shares the environment with the host. Whether we are using Docker or Singularity, the actions going on internally with the scientific filesystem client are the same. Given a simple enough pipeline, we could stop here, and just issue a series of commands to run the different apps.

### Deploy a container
Once you have your container completely functional you can upload it to Docker Hub or Singularity Hub for version control and sharing. We are uploading it to Docker Hub since Singularity can easily build or run a container from there.

If you have build your docker image in your container you can do:
```
# See your images
docker images
# You have to create an account on Docker Hub and a new repository which is very easy.
# Then log to docker hub from the command line.
docker login --username=yourhubusername --email=youremail@company.com
# You can tag your image.
docker tag bb38976d03cf buisciii/nextflow-scif:1.0
# You can push
docker push bu-isciii/nextflow-scif
```

## 2. Nextflow workflow
### Config file

### Main nextflow script

### Running the pipeline

## 3. Testing
**TODO**


# Advance usage
Till this point everything is working fine because we are using the data located in our current directory, and singularity mounts ```$PWD``` by default. But what if we are running in our HPC-cluster, and our data is in /processing_Data or other custom directory?

**TODO**: Try in different locations, with data and reference genome in different folders.
Use of --bind ? . Need to pre-create paths inside the containers??
Improvement of I/O metadata mounting only the project folder instead of all processing_Data?? 
*Note*: for each app scif/app/{app}/data folder is created. Thats a possibility for mounting. We have to make some performance test in our cluster.

# Summary
- We have created a container with all the software needed for our pipeline with Singularity and Docker using the same SCI-F recipes. This will allow you to deploy new containers really quicky for new pipelines once you have created a good collection of scif apps recipes.
- Using SCI-F we can access all our apps installed separately from our container. Besides our container is no longer a "black-box" and we can obtain information about what apps are installed in our container, with its labels and help.
- Our container is also functional with the native container commands, and we can exec any command in our container without the use of SCI-F which made it flexible enough for its use with nextflow.
- We have deploy our container to DockerHub so new users don't have to build the image theirselves.
- Finally we have created a basic nextflow workflow, with three different profile configurations {standard, singularity, docker}, this way you can use the same pipeline running in your local computer with your software installed, or with the singularity or docker image. Besides, you can create more profiles for its usage in a HPC with ```module``` configured, or even in a AWS server.

The idea is to have a workflow independent of where you are computing it, and a container system which allows you to assure reproducibility and sharing, including the maximun information possible for its use.

# Notes
For get this working I have to make a decision: put SCIF on the "top" or put nextflow instead. I will try to explain it:

- Nextflow on the top:
    * Nextflow provides a natively functionality for using singularity and docker. Using it you just set the commands the pipeline is going to run, and nextflow handles the use of ```docker run``` or ```singularity exec``` before this commands with the container you supplied in the configuration.
    * This means that the executables for each app **MUST** be accesible when you run/exec from the container, not each of the apps with SCI-F. 
    * Thus, we are being forced to add the executables to the PATH; and somehow mess up SCIF goal of separate apps environment adding the bin app folder to the PATH only when the app is run. Being this the main con of this approximation.
    * On the other hand, our main.nf will remain untouch and will be able to run in a local machine, a cluster or a container without any modification. Which was I was looking for.
    * For this to work your minimun needs in your system are nextflow and Docker/Singularity if you are using docker or singularity profile.
    
- SCIF on the top (I am going to make a branch testing this):
   * Another aproximation will be to install nextflow inside the containers. This way you can only interact with the apps commands through SCIF which will be the only entrypoint in Singularity and Docker.
   * Nextflow in this case will not use Singularity or Docker profile, because it is INSIDE one, so this config parameters must be set to false.
   * Main.nf script have to be modified and ```scif exec``` must be add before each command.
   * The pros of this approximation will be that apps will be in separate environments as SCI-F plans to, and you won't need to have nextflow installed in your system. If you use the image, just with docker or singularity installed you are set.
   * However you could only use this nextflow pipeline with the image prepare for this, or have scif installed in your local/hpc system, and all the apps you are going to use inside scif which may not be the case (as in our HPC where we don't have install permissions)
   
I suppose this is a decision you have to make, one could think that I can use my singularity image in my HPC system and forgot about modules and stuff, I will have all the software I need, in the version I need and whenever I need it. However even if this is going to be the case (we are still testing performance, we are very newbies using containers) I always like to have flexibility, and nextflow provides me in theory that "I don't mind where I compute" flexibility.
