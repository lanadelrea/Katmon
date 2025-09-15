# Start from micromamba base image
FROM mambaorg/micromamba:1.5.8

# Set working directory
WORKDIR /opt

# Copy environment YAML into container
COPY ~/Katmon/envs/jvarkit.yml .

# Create environment named "jvarkit"
RUN micromamba env create -f jvarkit.yml -y && \
    micromamba clean --all --yes

# Make sure the environment is activated by default
SHELL ["micromamba", "run", "-n", "jvarkit", "/bin/bash", "-c"]

# Add to path
ENV PATH="/opt/conda/envs/jvarkit/bin:$PATH"

# Default command: start bash inside the environment
CMD ["/bin/bash"]