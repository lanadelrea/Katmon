# Start with the specified micromamba base image
FROM mambaorg/micromamba:git-bb82854-lunar

# Set environment variables
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV MAMBA_USER=mambauser
ENV MAMBA_USER_ID=57439
ENV MAMBA_USER_GID=57439
ENV MAMBA_EXE=/bin/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV ENV_NAME=base
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Copy necessary scripts and setup
COPY _dockerfile_shell.sh /usr/local/bin/_dockerfile_shell.sh
COPY _activate_current_env.sh /usr/local/bin/_activate_current_env.sh
COPY _entrypoint.sh /usr/local/bin/_entrypoint.sh
COPY _apptainer_shell.sh /usr/local/bin/_apptainer_shell.sh
COPY _dockerfile_setup_root_prefix.sh /usr/local/bin/_dockerfile_setup_root_prefix.sh
COPY _dockerfile_initialize_user_accounts.sh /usr/local/bin/_dockerfile_initialize_user_accounts.sh
COPY /mnt/c/Pipelines/VirStrain/Custom_DB /Custom_DB

# Set shell to be used by Dockerfile
SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

# Set working directory
WORKDIR /tmp

# Set user
USER mambauser

# Install virstrain and clean up
RUN micromamba install -y -c bioconda -c conda-forge virstrain && \
    micromamba clean --all --yes

# Copy the new database file into the Docker image
COPY your_database_file.db /opt/virstrain/db/your_database_file.db

# Default command
CMD ["/bin/bash"]

# Entrypoint
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
