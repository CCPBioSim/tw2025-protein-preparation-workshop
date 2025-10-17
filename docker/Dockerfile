# Start with BioSim base image.
ARG BASE_IMAGE=latest
FROM ghcr.io/ccpbiosim/jupyterhub-base:$BASE_IMAGE

LABEL maintainer="James Gebbie-Rayet <james.gebbie@stfc.ac.uk>"

# Switch to jovyan user.
USER $NB_USER
WORKDIR $HOME

# Install workshop deps
RUN mamba install conda-forge::ambertools bioconda::blast conda-forge::openmm numpy=2.2 -y
RUN pip install git+https://github.com/CompBioAsia/CBAtools.git

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
