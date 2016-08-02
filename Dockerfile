FROM limesbonn/kallisto

MAINTAINER "Brett Beaulieu-Jones" brettbe@med.upenn.edu

# Libraries for processing + quantifying RNA-Seq
RUN cd /docker/
RUN apt-get install -y wget
RUN wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/curren$
RUN tar -vxzf sratoolkit.tar.gz
RUN export PATH=$PATH:$PWD/sratoolkit.2.4.0-1.mac64/bin

# Install helpful R libraries
RUN echo 'source("http://bioconductor.org/biocLite.R")' > /tmp/packages.R
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("rhdf5", "biomaRt", $
RUN echo 'install.packages("devtools")' > /tmp/packages.R
RUN echo 'devtools::install_github("pachterlab/sleuth")' > /tmp/packages.R
RUN echo 'library("sleuth")' > /tmp/packages.R
RUN echo 'install.packages("ggplot2")' > /tmp/packages.R
RUN echo 'library("ggplot2")' > /tmp/packages.R
RUN echo 'install.packages("plyr")' > /tmp/packages.R
RUN echo 'library("plyr")' > /tmp/packages.R
RUN echo 'install.packages("RColorBrewer")' > /tmp/packages.R
RUN echo 'library("RColorBrewer")' > /tmp/packages.R
RUN echo 'install.packages("stringr")' > /tmp/packages.R
RUN echo 'library("stringr")' > /tmp/packages.R
RUN echo 'install.packages("plyr")' > /tmp/packages.R
RUN echo 'library("plyr")' > /tmp/packages.R

RUN Rscript /tmp/packages.R

