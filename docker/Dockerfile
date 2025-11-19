# Start with BioSim base image.
ARG BASE_IMAGE=latest
FROM ghcr.io/ccpbiosim/jupyterhub-base:$BASE_IMAGE

LABEL maintainer="James Gebbie-Rayet <james.gebbie@stfc.ac.uk>"

ARG AMBERTOOLS_DL=null
ARG AMBERTOOLS_VER=25

# Install workshop deps
RUN mamba install -c conda-forge bioconda::blast openmm numpy=2.2 matplotlib scipy -y
RUN pip install git+https://github.com/CharlieLaughton/Alphafix.git
WORKDIR /tmp

RUN wget --no-verbose $AMBERTOOLS_DL/ambertools$AMBERTOOLS_VER.tar.bz2 && \
    tar xjf ambertools$AMBERTOOLS_VER.tar.bz2 && \
    rm ambertools$AMBERTOOLS_VER.tar.bz2 && \
    mv ambertools${AMBERTOOLS_VER}_src ambertools_src

WORKDIR /tmp/ambertools_src
RUN ./update_amber --update
RUN mkdir /tmp/build
WORKDIR /tmp/build

RUN cmake /tmp/ambertools_src -DCMAKE_INSTALL_PREFIX=/opt/conda -DPYTHON_EXECUTABLE=/opt/conda/bin/python -DUSE_CONDA_LIBS=TRUE -DBUILD_PYTHON=TRUE -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE -DDOWNLOAD_MINICONDA=FALSE -DCOMPILER=MANUAL -DBUILD_GUI=FALSE -DCOMPILER=GNU -DOPENMP=TRUE -DCUDA=FALSE -DMPI=FALSE -DUSE_FFT=True -DBUILD_DEPRECATED=False -DBUILD_INDEV=False -DBUILD_PERL=True -DOPTIMIZE=True
RUN make -j8
RUN make install

# Cleanup.
RUN rm -r /tmp/ambertools_src /tmp/build

# Switch to jovyan user.
USER $NB_USER
WORKDIR $HOME

# Add the exec to the path.
ENV AMBERHOME=/opt/conda

RUN mkdir blastp
WORKDIR /home/jovyan/blastp
RUN update_blastdb.pl --decompress swissprot
WORKDIR /home/jovyan

ENV BLASTP=/home/jovyan/blastp

# Get workshop files and move them to jovyan directory.
COPY --chown=1000:100 . .
RUN rm -rf INSTALL.txt LICENSE README.md docker .git .github

# Copy updated lab workspace
COPY --chown=1000:100 docker/default-37a8.jupyterlab-workspace /home/jovyan/.jupyter/lab/workspaces/default-37a8.jupyterlab-workspace

# Always finish with non-root user as a precaution.
USER $NB_USER
