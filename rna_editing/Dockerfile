#// Pull base image.
FROM centos

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/rna-editing
LABEL org.opencontainers.image.authors Yang Yang (yang.yang@bjmu.edu.cn)
LABEL org.opencontainers.image.title ICGC ARGO RNA-Seq editing image

#// PREPARE REDIportal environment: https://github.com/BioinfoUNIBA/REDItools/blob/master/Docker/Dockerfile
#Dockerfile Mantainer
LABEL mantainer="clalogiudice@gmail.com"

#Update the centos software with yum package-manager
RUN yum update -y && yum clean all

#Install git, wget and nano package
RUN yum -y install git wget nano && yum clean all

#Clone Nature_protocol Git repository
RUN git clone https://github.com/BioinfoUNIBA/REDItools

WORKDIR "/REDItools/NPscripts/" 

#Install miniconda with conda packages required by the nature_protocol
RUN chmod +x conda_pckg_installer_docker.py
RUN ./conda_pckg_installer_docker.py

#PREPARE NATURE_PROTOCOL environment
WORKDIR "/"
RUN echo "python ./REDItools/NPscripts/download-prepare-data-NP_docker.py" >> /root/.bashrc

#// install samtools, bcftools
WORKDIR "/opt"
RUN wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
RUN tar jxvf samtools-1.12.tar.bz2
RUN cd /opt/samtools-1.12/htslib-1.12
RUN ./configure --prefix=/opt/samtools-1.12/htslib-1.12/build
RUN make
RUN make install

WORKDIR "/opt/samtools-1.12"
RUN ./configure --prefix=/opt/samtools-1.12/build
RUN make
RUN make install
RUN ln -sf /opt/samtools-1.12/samtools ./

WORKDIR /opt/bcftools-1.12
RUN ./configure --prefix=/opt/bcftools-1.12/build --with-htslib=/opt/samtools-1.12/htslib-1.12
RUN ln -sf /opt/bcftools-1.12/bcftools ./

### set up environment
ENV LD_LIBRARY_PATH /opt/samtools-1.12/htslib-1.12/build/lib:$LD_LIBRARY_PATH
ENV PATH /miniconda2/bin:$PATH
#RUN echo "source activate nature_protocol" >> ~/.bashrc

ENV PATH="$PATH:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:"

CMD ["/bin/bash"]

