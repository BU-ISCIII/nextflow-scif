FROM centos

# Dependencies
RUN	yum -y groupinstall "Development Tools"
RUN yum -y update && yum -y install wget curl
RUN yum -y install python-setuptools
RUN easy_install pip

# Install scif from pypi
RUN pip install scif

# Install the filesystem from the recipe
ADD apps_recipes/*.scif /opt/
RUN scif install /opt/bwa_v0.7.17_centos7.scif
ENV PATH=${PATH}:/scif/apps/bwa/bin
RUN scif install /opt/samtools_v1.9_centos7.scif
ENV PATH=${PATH}:/scif/apps/samtools/bin
RUN scif install /opt/bcftools_v1.9_centos7.scif
ENV PATH=${PATH}:/scif/apps/bcftools/bin

# SciF Entrypoint
# Disabled because of compatibility with nextflow.
#ENTRYPOINT ["scif"]
#CMD scif

RUN find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\$PATH:$in" >> /etc/bashrc;done 
RUN if [ -z "${LD_LIBRARY_PATH-}" ]; then echo "export LD_LIBRARY_PATH=/usr/local/lib" >> /etc/bashrc;fi
RUN find /scif/apps -maxdepth 2 -name "lib" | while read in; do echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$in" >> /etc/bashrc;done
